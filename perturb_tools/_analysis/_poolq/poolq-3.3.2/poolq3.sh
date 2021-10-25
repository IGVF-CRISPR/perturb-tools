#!/bin/bash -ue

# This function figures out where a chain of symlinks points to:
function deref() {
    file="${1%/}"
    count=0
    while [ -L "$file" ]; do
        target=`readlink "$file"`
        [ ${target:0:1} == "/" ] || target=`dirname "$file"`/"$target"

        # strip trailing slashes; -L cannot handle those
        file="${target%/}"

        count=$(( count + 1 ))
        [ $count -eq 10 ] && exit 1
    done
    echo $file
}

programLoc="`dirname \`deref "$0"\``"

exec java -Xmx4G -XX:+UseG1GC -jar "$programLoc"/poolq3.jar ${1+"$@"}
