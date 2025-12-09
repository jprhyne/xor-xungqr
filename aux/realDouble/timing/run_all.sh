#!/bin/env bash

# List all the shell
# Source - https://stackoverflow.com/a
# Posted by Bennet Yee, modified by community. See post 'Timeline' for change history
# Retrieved 2025-12-08, License - CC BY-SA 4.0

find * -maxdepth 1 -type f -name "*.sh" -exec ./runFile {} $(basename "$0")\
 \;

echo "Finished!"
