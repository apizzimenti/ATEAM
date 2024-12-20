# Generate documentation with pdoc; rsync to the immediate parent directory;
# delete the old documentation directory.
pdoc ateam --force --html --template-dir docs/templates --output-dir=docs
rsync -a docs/ateam/ docs/
rm -rf docs/ateam


# Get the current path to this script: borrowed from https://bit.ly/3vHDuR4.
SCRIPT_HOME=`dirname $0 | while read a; do cd $a && pwd && break; done`

# Serve the documentation so we can take a look.
open "file://$SCRIPT_HOME/docs/index.html"