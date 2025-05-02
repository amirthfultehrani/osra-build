#!/bin/bash
mkdir -p $HOME/osra/2.2.1 || { echo "Cannot create $HOME/osra folder" 1>&2; exit; } 
cp --remove-destination package/* $HOME/osra/2.2.1/ || { echo "Cannot copy to $HOME/osra folder" 1>&2; exit; }
echo "Installing binary files in $HOME/osra"

