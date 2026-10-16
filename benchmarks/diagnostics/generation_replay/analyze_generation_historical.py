"""Independent target-level analysis of frozen campaigns; incomplete grids fail closed.

This evaluator does not import the project's aggregate_runs implementation.
Candidate-level similarity is conditional on successful scoring and is not activity.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import numpy as np

BASE = Path(__file__).resolve().parent
GROUPS = ['baselines', 'gpu0', 'gpu1']
SCORES = ['hit_shape', 'hit_pharmacophore', 'score']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(run):
    return (run['method_id'], run['query_id'], run.get('seed', run.get('generation_seed')))


def mean(values):
    vals = [v for v in values if v is not None]
    return statistics.mean(vals) if vals else None


def bootstrap(values):
    a = np.array([v for v in values if v is not None], dtype=float)
    result = {'mean': float(a.mean()) if len(a) else None, 'independent_targets': len(a), 'ci95': None}
    if len(a) > 1:
        draws = np.random.default_rng(20261016).choice(a, (10000, len(a)), replace=True).mean(axis=1)
        result['ci95'] = np.quantile(draws, [.025, .975]).tolist()
    return result


INFRASTRUCTURE_STATUS = 'infrastructure_interruption'
BOUNDED_METRICS = (
    'attempt_with_hop',
    'attempt_with_verified_route',
    'common_3d_output_slot_coverage',
)
# These describe the realized final outputs, not the unobserved chemistry outcome.
# The reference population and conditional scores must never be fabricated as zero.
EMPTY_OUTPUT_METRICS = (
    'run_success', 'candidate_count', 'valid_candidate_count', 'valid_fraction',
    'unique_families', 'duplicate_chemical_states', 'unchanged_query_families',
    'hop_families', 'unique_hop_scaffolds', 'unique_hop_ring_skeletons',
    'recovered_active_scaffolds', 'active_scaffold_precision_proxy',
    'exact_known_active_matches', 'known_inactive_family_matches',
    'decoy_family_matches', 'verified_rule_routes',
)


def availability(rows, targets):
    by_query = {q: sum(r['query_id'] == q for r in rows) for q in targets}
    by_target = {
        t: sum(by_query[q] for q in targets if targets[q] == t)
        for t in sorted(set(targets.values()))
    }
    return {
        'attempts': len(rows),
        'queries_with_observed_attempts': sum(n > 0 for n in by_query.values()),
        'targets_with_observed_attempts': sum(n > 0 for n in by_target.values()),
        'per_query_attempts': by_query,
        'per_target_attempts': by_target,
    }


def summarize_rows(rows, targets, names=None):
    """Preserve the original seed -> query -> target hierarchy, skipping only None."""
    names = sorted(names if names is not None else {n for r in rows for n in r['metrics']})
    per_query = {
        q: {m: mean(r['metrics'].get(m) for r in rows if r['query_id'] == q) for m in names}
        for q in targets
    }
    per_target = {
        target: {
            m: mean(per_query[q][m] for q in targets if targets[q] == target)
            for m in names
        }
        for target in sorted(set(targets.values()))
    }
    return {
        'per_query': per_query,
        'per_target': per_target,
        'target_macro': {m: bootstrap(t[m] for t in per_target.values()) for m in names},
    }


def summarize_bounds(rows, targets):
    """Identification bounds, with the full seed denominators; these are not CIs."""
    per_query = {
        q: {
            metric: {side: mean(r['metrics'][metric][side] for r in rows if r['query_id'] == q)
                     for side in ('lower', 'upper')}
            for metric in BOUNDED_METRICS
        }
        for q in targets
    }
    per_target = {
        target: {
            metric: {side: mean(per_query[q][metric][side] for q in targets if targets[q] == target)
                     for side in ('lower', 'upper')}
            for metric in BOUNDED_METRICS
        }
        for target in sorted(set(targets.values()))
    }
    return {
        'per_query': per_query,
        'per_target': per_target,
        'target_macro': {
            metric: {side: mean(t[metric][side] for t in per_target.values())
                     for side in ('lower', 'upper')}
            for metric in BOUNDED_METRICS
        },
    }


def infrastructure_sensitivity(per_method, targets):
    methods, complete_cases, paired, paired_complete = {}, {}, {}, {}
    indexed = {}
    for method, rows in sorted(per_method.items()):
        names = sorted({m for r in rows for m in r['metrics']})
        observed = [r for r in rows if r['status'] != INFRASTRUCTURE_STATUS]
        bounded = []
        for row in rows:
            missing = row['status'] == INFRASTRUCTURE_STATUS
            values = {}
            for metric in BOUNDED_METRICS:
                value = row['metrics'][metric]
                if not isinstance(value, (int, float)) or not np.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f'{method}/{row["query_id"]}: {metric} outside its [0,1] domain')
                values[metric] = {'lower': 0.0 if missing else value, 'upper': 1.0 if missing else value}
            bounded.append({**row, 'metrics': values})
        indexed[method] = {(r['query_id'], r['seed']): r for r in bounded}
        if len(indexed[method]) != len(rows):
            raise ValueError('Duplicate query/seed within method')
        methods[method] = {
            'attempts': len(rows),
            'infrastructure_missing_attempts': len(rows) - len(observed),
            'observed_outcome_attempts': len(observed),
            'full_grid_availability': availability(rows, targets),
            'observed_outcome_availability': availability(observed, targets),
            **summarize_bounds(bounded, targets),
        }
        complete_cases[method] = {
            'attempts': len(observed),
            'excluded_infrastructure_attempts': len(rows) - len(observed),
            'status_counts': dict(Counter(r['status'] for r in observed)),
            'availability': availability(observed, targets),
            'runtime_observations': availability(
                [r for r in observed if r['metrics']['wall_seconds'] is not None], targets),
            **summarize_rows(observed, targets, names),
        }
    if 'pipeline' in indexed:
        primary = indexed['pipeline']
        for method, comparison in indexed.items():
            if method == 'pipeline':
                continue
            if primary.keys() != comparison.keys():
                raise ValueError('Paired sensitivity requires the same attempted query/seed grid')
            bounds, complete = [], []
            both_missing = any_missing = 0
            for pair_key, left in primary.items():
                right = comparison[pair_key]
                left_missing = left['status'] == INFRASTRUCTURE_STATUS
                right_missing = right['status'] == INFRASTRUCTURE_STATUS
                any_missing += left_missing or right_missing
                both_missing += left_missing and right_missing
                difference = {
                    metric: {
                        'lower': left['metrics'][metric]['lower'] - right['metrics'][metric]['upper'],
                        'upper': left['metrics'][metric]['upper'] - right['metrics'][metric]['lower'],
                    }
                    for metric in BOUNDED_METRICS
                }
                bounds.append({**left, 'metrics': difference})
                if not left_missing and not right_missing:
                    complete.append({**left, 'metrics': {m: value['lower'] for m, value in difference.items()}})
            paired[method] = {
                'paired_attempts': len(bounds),
                'pairs_with_infrastructure_missing_outcome': any_missing,
                'pairs_with_both_outcomes_missing': both_missing,
                **summarize_bounds(bounds, targets),
            }
            paired_complete[method] = {
                'paired_attempts': len(complete),
                'excluded_infrastructure_pairs': any_missing,
                'availability': availability(complete, targets),
                **summarize_rows(complete, targets, BOUNDED_METRICS),
            }
    return {
        'infrastructure_missing_outcome_sensitivity': {
            'interpretation': 'Infrastructure interruptions are missing counterfactual outcomes, not chemistry failures. Bounds range each missing bounded attempt over [0,1] while retaining all planned seed/query/target denominators. They are identification bounds, not confidence intervals or imputed results.',
            'bounded_metric_domains': {m: [0, 1] for m in BOUNDED_METRICS},
            'hierarchy': 'Average all attempted seeds within each query, queries within each target, then equal target weights. Paired bounds subtract comparator upper from pipeline lower, and comparator lower from pipeline upper at matching query/seed slots.',
            'methods': methods,
            'paired_pipeline_minus_comparator': paired,
            'other_metrics': 'No missing-outcome upper bound is claimed for unbounded metrics such as runtime. Other metrics remain realized-output or success-conditional descriptions; no counterfactual upper bounds are invented for them.',
        },
        'complete_case_secondary': {
            'interpretation': 'Secondary, selection-biased availability analysis: exclude only infrastructure_interruption. Timeouts and other recorded method failures remain zero-output attempts. Method means use available seeds/queries/targets; paired differences require both outcomes at the same query/seed. Availability and denominators can differ, so these are not estimates under random missingness.',
            'methods': complete_cases,
            'paired_pipeline_minus_comparator_bounded_metrics': paired_complete,
        },
    }


def analyze_runs(runs, scored, targets, caps):
    """Pure aggregation entry point; synthetic tests do not touch campaign outputs."""
    score_lookup = {(key(s), s['rank'], s['smiles']): s for s in scored}
    if len(score_lookup) != len(scored):
        raise ValueError('Duplicate common-score occurrence')
    reference_counts = {}
    for query in targets:
        counts = {
            r.get('evaluation', {}).get('metrics', {}).get('reference_active_scaffolds')
            for r in runs if r['query_id'] == query and r['status'] != INFRASTRUCTURE_STATUS
        } - {None}
        if len(counts) > 1:
            raise ValueError('Inconsistent reference population for a query')
        reference_counts[query] = next(iter(counts), None)
    per_method = defaultdict(list)
    occurrence_count = 0
    for run in runs:
        if run['target_id'] != targets[run['query_id']]:
            raise ValueError('Run target does not match query target')
        infrastructure = run['status'] == INFRASTRUCTURE_STATUS
        evaluation = run.get('evaluation', {})
        candidates = evaluation.get('candidates', []) if infrastructure else evaluation['candidates']
        if infrastructure and candidates:
            raise ValueError('Infrastructure missing-outcome marker must not contain accepted candidates')
        occurrence_count += len(candidates)
        rows = [score_lookup[(key(run), rank, c['smiles'])] for rank, c in enumerate(candidates, 1)]
        # Do not turn missing similarity into a zero score.
        success = [r for r in rows if r['all_3d_seeds_succeeded']]
        if infrastructure:
            metrics = {m: 0 for m in EMPTY_OUTPUT_METRICS}
            reference_count = reference_counts[run['query_id']]
            metrics['reference_active_scaffolds'] = reference_count
            metrics['active_scaffold_recall_proxy'] = 0.0 if reference_count else None
        else:
            metrics = dict(evaluation['metrics'])
        runtime = run['resources']['wall_seconds']
        if runtime is not None and (not isinstance(runtime, (int, float)) or not np.isfinite(runtime) or runtime < 0):
            raise ValueError('Runtime must be a finite nonnegative observation or null')
        metrics.update({
            'attempt_with_hop': int(metrics['hop_families'] > 0),
            'attempt_with_verified_route': int(metrics['verified_rule_routes'] > 0),
            'common_3d_successful_candidates': len(success),
            'common_3d_failed_candidates': len(rows) - len(success),
            'common_3d_output_slot_coverage': len(success) / caps[run['method_id']],
            'wall_seconds': runtime,
        })
        unique = {}
        for candidate in candidates:
            unique.setdefault(candidate['family'], candidate)
        support = [c for c in unique.values() if c.get('declared_external_synthesis_evidence')]
        metrics['declared_forward_assessed_families'] = sum(
            bool(c['declared_external_synthesis_evidence'].get('steps')) and
            all(s['forward_prediction']['status'] in ['matched', 'mismatched']
                for s in c['declared_external_synthesis_evidence']['steps']) for c in support)
        metrics['declared_forward_all_steps_matched_families'] = sum(
            bool(c['declared_external_synthesis_evidence'].get('steps')) and
            all(s['forward_prediction']['status'] == 'matched'
                for s in c['declared_external_synthesis_evidence']['steps']) for c in support)
        for metric in SCORES:
            metrics[f'common_{metric}_mean_success_conditional'] = mean(r['scores'][metric] for r in success)
            metrics[f'common_{metric}_range_mean_success_conditional'] = mean(
                r['scores'][metric + '_stability']['range'] for r in success)
        per_method[run['method_id']].append({
            'query_id': run['query_id'], 'target_id': run['target_id'],
            'seed': key(run)[2], 'status': run['status'], 'metrics': metrics,
        })
    if occurrence_count != len(scored):
        raise ValueError('Missing or extraneous scored occurrence')
    results = {}
    for method, rows in sorted(per_method.items()):
        runtime_rows = [r for r in rows if r['metrics']['wall_seconds'] is not None]
        results[method] = {
            'attempts': len(rows),
            'status_counts': dict(Counter(r['status'] for r in rows)),
            'infrastructure_interruption_attempts': sum(r['status'] == INFRASTRUCTURE_STATUS for r in rows),
            'runtime_observed_attempts': len(runtime_rows),
            'runtime_missing_attempts': len(rows) - len(runtime_rows),
            'runtime_measurement_status_counts': dict(Counter(
                r['resources'].get('measurement_status', 'observed' if r['resources']['wall_seconds'] is not None else 'missing_unspecified')
                for r in runs if r['method_id'] == method)),
            'runtime_observation_availability': availability(runtime_rows, targets),
            **summarize_rows(rows, targets),
        }
    differences = {}
    if 'pipeline' in results:
        primary = results['pipeline']['per_target']
        for method in results:
            if method == 'pipeline':
                continue
            comparison = results[method]['per_target']
            differences[method] = {}
            for metric in results['pipeline']['target_macro']:
                vals = [primary[t][metric] - comparison[t].get(metric) for t in primary
                        if primary[t][metric] is not None and comparison[t].get(metric) is not None]
                differences[method][metric] = bootstrap(vals)
    return {
        'methods': results,
        'paired_pipeline_minus_comparator': differences,
        **infrastructure_sensitivity(per_method, targets),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--split', choices=['scaffold', 'time', 'target'], required=True)
    parser.add_argument('--allow-partial', action='store_true', help='Write grid audit only; do not compare unfinished methods')
    args = parser.parse_args()
    split = args.split
    panel_path = BASE / f'panels/direct-human-{split}/panel.json'
    panel = json.loads(panel_path.read_text())
    records = {r['id']: r for r in panel['compounds']}
    targets = {q['id']: records[q['compound_id']]['target_id'] for q in panel['queries']}
    expected, runs, sources, caps = set(), [], {}, {}
    for group in GROUPS:
        pp = BASE / (f'main-protocol-{group}.json' if split == 'scaffold' else f'secondary-{split}-protocol-{group}.json')
        protocol = json.loads(pp.read_text())
        sources[str(pp.relative_to(BASE))] = digest(pp)
        for method in protocol['methods']:
            mid = method['id']
            caps[mid] = method.get('budget_overrides', {}).get('max_candidates', protocol['budget']['max_candidates'])
            expected.update((mid, q, seed) for q in targets for seed in protocol['seeds'])
        folder = f'comparison-{group}' if split == 'scaffold' else f'comparison-{split}-{group}'
        rp = BASE / folder / 'runs.json'
        if rp.exists():
            runs.extend(json.loads(rp.read_text()))
            sources[str(rp.relative_to(BASE))] = digest(rp)
    actual = Counter(key(r) for r in runs)
    missing = sorted(expected - actual.keys())
    unexpected = sorted(actual.keys() - expected)
    duplicates = [list(k) for k, n in actual.items() if n > 1]
    for r in runs:
        assert r['target_id'] == targets[r['query_id']]
    out = BASE / f'analysis-generation-{split}'
    out.mkdir(exist_ok=True)
    audit = {'split': split, 'expected_runs': len(expected), 'observed_runs': len(runs),
             'missing': missing, 'unexpected': unexpected, 'duplicates': duplicates,
             'status_counts': dict(Counter(r['status'] for r in runs)),
             'complete': not (missing or unexpected or duplicates)}
    (out / 'grid-audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    if not audit['complete']:
        print(json.dumps({k: v for k, v in audit.items() if k not in ['missing', 'unexpected', 'duplicates']}))
        if args.allow_partial:
            return
        raise SystemExit('Incomplete/invalid attempt grid; no comparative estimates emitted')
    score_dir = 'common-generated-scores-8c' if split == 'scaffold' else f'common-generated-{split}-scores-8c'
    sp = BASE / score_dir / 'generated_scores.json'
    progress = json.loads((sp.parent / 'progress.json').read_text())
    if not progress['complete']:
        raise SystemExit('Independent common scoring incomplete')
    scored = json.loads(sp.read_text())
    sources[str(sp.relative_to(BASE))] = digest(sp)
    analysis = analyze_runs(runs, scored, targets, caps)
    report = {
        'split': split, 'grid_audit': audit, 'source_sha256': sources, 'dataset_sha256': digest(panel_path),
        **analysis,
        'uncertainty': 'Average seeds within query, queries within target, then equal targets. 10,000 target bootstrap draws, seed 20261016. Descriptive unadjusted intervals; no confirmatory significance or multiplicity claim.',
        'scope': 'All attempted inputs retained, including timeouts, zero-output runs and explicit infrastructure interruptions. Primary methods/differences describe conservative realized final output; an infrastructure missing outcome is not classified as chemistry failure. Separate missing-outcome bounds and selection-biased complete-case analyses accompany that view. Common similarity means remain conditional on score success; null runtime is excluded from runtime means with observed counts reported. Declared forward-model support is not independent re-inference. Structural recovery is a proxy, not predicted experimental activity or yield.',
        'assessment_design_note': 'Independent aggregation specified during the ongoing campaign; generation methods and budgets remain frozen. Repeated targets across split panels are not new independent targets.'}
    (out / 'comparison.json').write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(json.dumps({'split': split, 'attempts': len(runs), 'methods': len(analysis['methods']), 'targets': len(set(targets.values())), 'complete': True}))


if __name__ == '__main__':
    main()
