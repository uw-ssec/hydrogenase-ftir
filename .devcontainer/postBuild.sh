# For writing commands that will be executed after the container is created

# Install the Bruker Opus Reader
python3 -m pip install matplotlib brukeropusreader openpyxl
python3 -m pip install -e .
mkdir -p "/workspaces/hydrogenase-ftir/data"