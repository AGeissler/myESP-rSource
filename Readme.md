# ESP-r: A dynamic building performance simulation program.

## 6. INSTALLATION

MacOS Monterey 12.6 :

XCode:
14.0.1
October 1st, 2022

achim@MacBook-Pro BT_Augustin % pkgutil --pkg-info=com.apple.pkg.CLTools_Executables
package-id: com.apple.pkg.CLTools_Executables
version: 14.0.0.0.1.1661618636
volume: /
location: /
install-time: 1663219055
groups: com.apple.FindSystemFiles.pkg-group 

Using gcc-md-7

Add 
(export) LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"

From:
https://stackoverflow.com/questions/56156520/gfortran-error-ld-library-not-found-for-lsystem-when-trying-to-compile 

Search in google:
  macos gcc-md ld: library not found for "-lSystem”

