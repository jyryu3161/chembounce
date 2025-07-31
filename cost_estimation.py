"""
Cost estimation: time, memory, CPU
"""
import datetime
    
import time
import psutil
import threading
from functools import wraps


class EstDeco:
    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __call__(self, *args, **kwargs):
        # Record start time and initial CPU times
        start_time = time.time()
        process = psutil.Process()
        cpu_times_before = process.cpu_times()
        system_cpu_times_before = psutil.cpu_times()

        # Execute the function
        result = self.func(*args, **kwargs)

        # Record end time and final CPU times
        end_time = time.time()
        cpu_times_after = process.cpu_times()
        system_cpu_times_after = psutil.cpu_times()

        # Calculate elapsed time
        elapsed_time = end_time - start_time

        # Calculate process CPU times
        user_cpu_time = cpu_times_after.user - cpu_times_before.user
        system_cpu_time = cpu_times_after.system - cpu_times_before.system
        total_cpu_time = user_cpu_time + system_cpu_time

        # Calculate system CPU times
        system_total_time_before = sum(system_cpu_times_before)
        system_total_time_after = sum(system_cpu_times_after)
        system_total_time = system_total_time_after - system_total_time_before

        # Calculate CPU usage percentage
        cpu_usage_percent = (total_cpu_time / system_total_time) * 100 * psutil.cpu_count()

        # Compile metrics into a dictionary
        cost = {
            'elapsed_time': elapsed_time,
            'user_cpu_time': user_cpu_time,
            'system_cpu_time': system_cpu_time,
            'total_cpu_time': total_cpu_time,
            'cpu_usage_percent': cpu_usage_percent,
            'cpu_count': psutil.cpu_count()
        }

        return result, cost
    

class EstDecoMem:
    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __call__(self, *args, **kwargs):
        # Record start time and initial CPU times
        start_time = time.time()
        process = psutil.Process()
        cpu_times_before = process.cpu_times()
        system_cpu_times_before = psutil.cpu_times()
        
        # Memory tracking variables
        memory_samples = []
        memory_monitoring = True
        
        # Get initial memory usage
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_samples.append(initial_memory)
        
        # Memory monitoring thread
        def monitor_memory():
            while memory_monitoring:
                try:
                    current_memory = process.memory_info().rss / 1024 / 1024  # MB
                    memory_samples.append(current_memory)
                    time.sleep(0.01)  # Sample every 10ms
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break
        
        # Start memory monitoring in a separate thread
        memory_thread = threading.Thread(target=monitor_memory, daemon=True)
        memory_thread.start()

        try:
            # Execute the function
            result = self.func(*args, **kwargs)
        finally:
            # Stop memory monitoring
            memory_monitoring = False
            memory_thread.join(timeout=0.1)

        # Record end time and final CPU times
        end_time = time.time()
        cpu_times_after = process.cpu_times()
        system_cpu_times_after = psutil.cpu_times()

        # Calculate elapsed time
        elapsed_time = end_time - start_time

        # Calculate process CPU times
        user_cpu_time = cpu_times_after.user - cpu_times_before.user
        system_cpu_time = cpu_times_after.system - cpu_times_before.system
        total_cpu_time = user_cpu_time + system_cpu_time

        # Calculate system CPU times
        system_total_time_before = sum(system_cpu_times_before)
        system_total_time_after = sum(system_cpu_times_after)
        system_total_time = system_total_time_after - system_total_time_before

        # Calculate CPU usage percentage
        cpu_usage_percent = (total_cpu_time / system_total_time) * 100 * psutil.cpu_count() if system_total_time > 0 else 0

        # Calculate memory statistics
        max_memory = max(memory_samples) if memory_samples else initial_memory
        avg_memory = sum(memory_samples) / len(memory_samples) if memory_samples else initial_memory
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_delta = final_memory - initial_memory

        # Compile metrics into a dictionary
        cost = {
            'elapsed_time': elapsed_time,
            'user_cpu_time': user_cpu_time,
            'system_cpu_time': system_cpu_time,
            'total_cpu_time': total_cpu_time,
            'cpu_usage_percent': cpu_usage_percent,
            'cpu_count': psutil.cpu_count(),
            'initial_memory_mb': initial_memory,
            'final_memory_mb': final_memory,
            'max_memory_mb': max_memory,
            'avg_memory_mb': avg_memory,
            'memory_delta_mb': memory_delta,
            'memory_samples_count': len(memory_samples)
        }

        return result, cost
    