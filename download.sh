wget https://zenodo.org/record/7810285/files/data.zip
unzip data.zip -d ./data
mv ./data/data/* ./data/
rm -r ./data/data/
rm data.zip