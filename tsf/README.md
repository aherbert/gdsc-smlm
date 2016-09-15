Support for TSF (Tagged Spot File) format
=========================================

The goal of the tagged spot file format is to provide an efficient data format 
for superresolution microscopy data that generate images by locating the 
position of single fluorescent emitters.

The orignal TSF description is available here:
[Tagged Spot File (TSF) format](https://micro-manager.org/wiki/Tagged_Spot_File_(tsf)_format)

The format has been extended to add support for fields in the GDSC SMLM
software.


Build
-----

The Google .proto buffer file is used to generate Java code using 
[Google Protocol buffers](https://developers.google.com/protocol-buffers/).

You must have the binary protoc executable in your path. You can generate
the Java code for teh GDSC SMLM project using the following command:

        protoc --java_out=../src/main/java TSFProto.proto

This will produce gdsc.smlm.tsf.TaggedSpotFile.java in the source tree.

This file can be opened and the following added before the TaggedSpotFile class 
to disable compiler warnings:

        @SuppressWarnings({"unchecked", "unused"})

Other edits can be performed at risk of breaking the code!
