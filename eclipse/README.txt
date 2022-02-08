The GDSC-SMLM code was developed using the Eclipse IDE:
https://eclipse.org/

The code can be built using Maven. See the README.md for details. However using
Eclipse is preferred for code development to provide a debugging environment.

You will need the Maven and Git Eclipse plugins. The standard Eclipse IDE for
Java developers has these.

Set-up the project
------------------

Import the project into Eclipse (File>Import)
Select: Maven > Existing Maven projects

This will import the project but may not link it to the source Git repository.
Right-click on the project name and select 'Team > Share'. If you share it back to
the same location it will attach to the source Git repository.

Code formatting
---------------

The Eclipse code format rules (in this directory) can be loaded using:

Eclipse Preferences : Java > Code Style > Formatter
Eclipse Preferences : Java > Code Style > Clean Up

Click 'Import...' to load the provided rules.

Running the code
----------------

Build the project.

---
Optional:

To debug the ImageJ 3D Viewer from Eclipse relies on being able to find the native
runtime libraries for Open GL. These should now be included by Maven. If this fails
then do a manual install:

1. Download JOGL
   http://jogamp.org/wiki/index.php/Downloading_and_installing_JOGL
2. Extract all the files
3. Copy jogl-all.jar and gluegen-rt.jar to a folder
4. Copy the 'natives' jars for the current platform to the same location,
   e.g. jogl-all-natives-linux-amd64.jar and gluegen-rt-natives-linux-amd64.jar
5. Add the jogl-all.jar and gluegen-rt.jar as external jars to the build path in the
   Eclipse project
---

Create a symbolic link on the filesystem to set-up the folders that are expected by ImageJ.

Windows:

    GDSC-SMLM>mklink /D gdsc-smlm-ij\plugins gdsc-smlm-ij\target\classes
    symbolic link created for plugins <<===>> target\classes

Linux:

    [GDSC-SMLM] % ln -s gdsc-smlm-ij/target/classes gdsc-smlm-ij/plugins

Create a new Run configuration in the gdsc-smlm-ij_ project.

Select ij.ImageJ as the main class.

Run the code.

The same target can be used to debug the code.
