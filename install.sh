wget https://zenodo.org/records/16741967/files/data.zip

unzip data.zip -d ./data
mv ./data/data/* ./data/
rm -r ./data/data/
rm data.zip
