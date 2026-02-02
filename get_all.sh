#!/bin/bash
set -x #echo on
remote_url=$(git config --get remote.origin.url)
for branch in $(git branch --all | grep '^\s*remotes' | egrep --invert-match '(:?HEAD|master)$'); do
		branch_name=$(echo $branch| cut -d'/' -f 3)
    git clone -b $branch_name $remote_url $branch_name
done

