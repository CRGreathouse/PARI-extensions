#!/usr/bin/sh

showChanges() {
	if [ "$same" -eq 1 ]; then
		# There are no changes in the CHANGES file
		echo "Nothing documented in CHANGES"
	else
		# There are changes in the CHANGES file
		echo "Documented changes:"
		git diff $oldHEAD $newHEAD -- CHANGES
	fi
	echo "Number of commits: $commits"
}

# Check tune.h
if [ -f Olinux-x86_64/tune.h ]; then
	lastTuned=$(date -r Olinux-x86_64/tune.h +"%b %d %H:%M")
	echo "Last tuned $lastTuned"
else
	echo "*** Warning: No tuning file detected; consider running"
	echo "./Configure --tune"
fi

# Pull the latest source from git, recording the respective HEADs.
oldHEAD=$(git rev-parse HEAD)
git pull
err=$?
if [ "$err" -ne 0 ]; then
	echo "Error $err with git pull, exiting."
	exit $err
fi
newHEAD=$(git rev-parse HEAD)

# If there are no updates, exit.
if [ "$newHEAD" = "$oldHEAD" ]; then
	last=$(git show -s --pretty=format:"%ah, %ar" $newHEAD)
	echo "  [no changes since $last]"
	exit 0
fi

# Show changes and number of commits
same=$(git rev-parse $oldHEAD:CHANGES $newHEAD:CHANGES | uniq | wc -l)
commits=$(git rev-list $newHEAD ^$oldHEAD --pretty=oneline --count)
showChanges

# Try to make, catching errors
make all
err=$?
if [ "$err" -ne 0 ]; then
	echo "Error $err with make all, exiting."
	echo "Possible fix: ./Configure --tune"
	exit $err
fi

# Check that the build was successful and, if so, install.
make dobench
err=$?
if [ "$err" -ne 0 ]; then
	echo "Error $err when running bench, build was unsuccessful!"
	exit $err
fi

sudo true
err=$?
if [ "$err" -ne 0 ]; then
	echo "Error $err with sudo; not able to install. Local build seems ok."
	echo "Check the error and consider retrying sudo make install."
	exit $err
fi

sudo make install
err=$?
if [ "$err" -ne 0 ]; then
	echo "Error $err when installing, but local build seems ok."
	echo "Check the error and consider retrying sudo make install."
	exit $err
fi

# Show the user the changes and number of commits again.
echo "Updated successfully. Documented changes:"
git diff $oldHEAD $newHEAD -- CHANGES
echo "Number of commits: $commits"
