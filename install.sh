#!/bin/bash
set -eu

pwd; hostname; date

echo "installing minizinc ..."
# Install MiniZinc
# Ask user for MiniZinc installation path
read -p "Enter the installation path for MiniZinc: " minizinc_path

# Append MiniZinc directory to the entered path
minizinc_path_full="$minizinc_path/MiniZinc"

if [ -d $minizinc_path_full ]; then
    rm -r $minizinc_path_full
fi

mkdir $minizinc_path_full
LATEST_MINIZINC_VERSION=$(curl -s https://api.github.com/repos/MiniZinc/MiniZincIDE/releases/latest | grep -oP '"tag_name": "\K(.*)(?=")')
wget "https://github.com/MiniZinc/MiniZincIDE/releases/download/$LATEST_MINIZINC_VERSION/MiniZincIDE-$LATEST_MINIZINC_VERSION-bundle-linux-x86_64.tgz"
tar -xvzf "MiniZincIDE-$LATEST_MINIZINC_VERSION-bundle-linux-x86_64.tgz" -C "$minizinc_path_full" --strip-components=1
rm "MiniZincIDE-$LATEST_MINIZINC_VERSION-bundle-linux-x86_64.tgz"
if [ -L /usr/local/bin/minizinc ]; then
    rm /usr/local/bin/minizinc
fi
ln -s "$minizinc_path_full/bin/minizinc" /usr/local/bin/minizinc

# Install Python requirements
python3 -m pip install -r pyrequirements.txt

date
echo "installation was completed!"
exit 0
