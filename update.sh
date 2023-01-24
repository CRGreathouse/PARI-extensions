#!/usr/bin/sh

# Pull the latest source from git, recording the respective HEADs.
oldHEAD=$(git rev-parse HEAD)
git pull
if [ $? -ne 0 ]; then
	echo "Error $? with git pull, exiting."
	exit 1
fi
newHEAD=$(git rev-parse HEAD)

# If there are no updates, exit.
if [ "$newHEAD" = "$oldHEAD" ]; then
	last=$(git show -s --pretty=format:"%ah, %ar" $newHEAD)
	echo "  [no changes since $last]"
	exit 0
fi

# Show changes and number of commits
echo "Documented changes:"
git diff $oldHEAD $newHEAD -- CHANGES
commits=$(git rev-list $newHEAD ^$oldHEAD --pretty=oneline --count)
echo "Number of commits: $commits"

# Try to make, catching errors
make all
if [ $? -ne 0 ]; then
	echo "Error $? with make all, exiting."
	echo "Possible fix: ./Configure --tune"
	exit 1
fi

# Check that the build was successful and, if so, install.
make dobench && sudo make install

# Show the user the changes and number of commits again.
echo "Updated successfully. Documented changes:"
git diff $oldHEAD $newHEAD -- CHANGES
echo "Number of commits: $commits"
