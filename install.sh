#!/bin/bash

# Try downloading from the first URL
echo "Attempting to download data from the first URL..."
if wget -O data.zip "https://www.dropbox.com/scl/fi/1wlp71fdvjycee8r52wp6/data.zip?rlkey=o328bgjyj2mtyf71khmzceh72&st=89du7ywu&dl=1"; then
    echo "Download successful from the first URL!"
else
    echo "First URL download failed. Trying the second URL..."
    if wget -O data.zip https://zenodo.org/records/16741967/files/data.zip; then
        echo "Download successful from the second URL!"
    else
        echo "Download failed from both URLs!"
        exit 1
    fi
fi

unzip data.zip -d ./data
mv ./data/data/* ./data/
rm -r ./data/data/
rm data.zip
