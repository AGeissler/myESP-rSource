#!/bin/bash
#
# This script copies new or changed files from a new official ESP-r
# release archive to the existing (current) source tree and deletes files
# found in the current release but missing in new release.
#
# Basic string manipulation from https://spin.atomicobject.com/2014/02/16/bash-string-maniuplation/
#
# 0. Step : Prepare the environment
IFS=$'\n'

# 1. Step : Set up some variables
TEST=false
#prefix="/Users/achim.geissler/Downloads"
prefix="/Users/achim/Downloads"
release="ESP-r_V13.3.9_Src"

# 2. Step : create file listings for the source root directory and
#           subdirectories src, models, modish, data and doc:
diff -q --exclude=CompileQuiet --exclude=diff_root.txt ./ ${prefix}/${release}/. > diff_root.txt
diff -rq src ${prefix}/${release}/src > diff_src.txt
diff -rq models ${prefix}/${release}/models > diff_models.txt
diff -rq modish ${prefix}/${release}/modish > diff_modish.txt
diff -rq data ${prefix}/${release}/data > diff_data.txt
diff -rq doc ${prefix}/${release}/doc > diff_doc.txt
#
#    The diffs are issued in the git repository directory "myESP-rSource"
#    (where this script should actually reside).
#
# 3. Step : clean the resulting diff files.
#    Following entry types are found:                          Following actions are desired:
#    "Only in <curpath>: <filename>"                           => delete from curpath
#    "Only in <newpath>: <filename>"                           => copy from newpath
#    "Files <curpath/filename> and <newpath/filename> differ"  => copy from newpath
#
#    Where <curpath> = <folder>[/<folder>[/...]]
#          <newpath> = <prefix>/<release>/<curpath>
#

filelist=(diff_src.txt diff_models.txt diff_modish.txt diff_data.txt diff_doc.txt diff_root.txt)
#filelist=(diff_models.txt)

for file in "${filelist[@]}"
do
   lines=`cat ${file}`
   for line in $lines
   do

#... if line begins with "Only" then
     if [[ $line == O* ]] ; then
#      remove "Only in "
       s1=${line#"Only in "}

#      split at ": " and then remove leading whitespace
       fn=${s1##*:}
       fn="${fn#"${fn%%[![:space:]]*}"}"

       pth=${s1%:*}
       pth2=${pth#"${prefix}/${release}/"*}

       len_pth=${#pth}
       len_pth2=${#pth2}

#      if part[1] contains <prefix>/<release> then
       if [ "$len_pth" -gt "$len_pth2" ] ; then

#        copy <prefix>/<release>/part[1]/part[2] to part[1]/part[2]
         if [ "$TEST" = "true" ]; then
           echo "cp -R ${pth}/${fn} ${pth2}/${fn}"
         else
           cp -R ${pth}/${fn} ${pth2}/${fn}
         fi
       else
#         if [ ${fn} != ".git" ]; then
         if [ ${fn} != ".git" ] && [ ${fn} != ".gitignore" ] \
                                && [ ${fn} != ".gitmodules" ] \
                                && [ ${fn} != ".project" ] \
                                && [ ${fn} != "release_copyfiles.sh" ]; then
           if [ "$TEST" = "true" ]; then
             echo "rm -R ${pth2}/${fn}"
           else
             rm -R ${pth2}/${fn}
           fi
         fi
       fi

#... else (line does not begin with "Only")
     elif [[ $line == F* ]] ; then
#      remove "Files "
       s1=${line#"Files "}
#      remove " differ"
       s2=${s1%" differ"}
#      split at " and "
       from=${s2#*" and "}
       to=${s2%" and "*}

#      copy part[2] to part[1]
       if [ "$TEST" = "true" ]; then
         echo "cp -R ${from} ${to}"
       else
         cp -R ${from} ${to}
       fi
     fi
#
   done # lines
done # filelist
