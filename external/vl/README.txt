Contents

    The VL vector & matrix package, version 1.3.2

    VL provides 2-, 3- and 4-vector and matrix types, arbitrarily-sized
    vector and matrix types, sparse and sub-region types, and lots of
    useful functions and arithmetic operators defined on those types. It
    includes a debugging version of the library, with range-checking of
    all indices.

    VL is free for commercial and non-commercial use; see the LICENSE
    file for redistribution conditions. (These apply only to the source
    code; binaries may be freely redistributed, no strings attached.)

Author

    Andrew Willmott, ajw+vl@cs.cmu.edu
    Please send me any suggestions/comments/bug fixes you have.

History

    Briefly: VL, and its simpler cousin SVL, have their roots in an
    attempt replace the various (inevitably heavily-macroized)
    vector/matrix packages used by the members of the CMU graphics group
    with a single C++ one, so that people could share code more easily.
    VL makes heavy use of inlines to achieve performance equal to that
    of the C macros and functions it replaced.
    
Compiling VL under Windows

    The distribution should come with workspace and project files for
    creating a static library, vl.lib. If these fail to work for some
    reason, below are the instructions used to create them.
    
    (These instructions are for MS Visual C++. Hopefully you can adapt
    them if you have a different windows compiler.)

    (a) Create a static library project called "vl" in the vl directory.

    (b) Add the files Basics.cxx, LibVLd.cxx, LibVLf.cxx and LibVLfd.cxx
    to the project via Project/Add To Project/Files.

    (c) In Project/Project Settings, add the "include" directory to the
    C/C++ pane under Category: Preprocessor, "Additional include
    directories".

    (d) If you want to create a debugging version of the library, add
    the symbol VL_CHECKING to the C/C++ pane's "Preprocesser
    definitions" field.

    Once you have the project open and set up correctly, hit F7 to build
    the library.

To Use

    (a) Add the vl.lib file created in the above section to your own
    project. You can use Project/Add To Project/Files to do this, though
    you'll have to change the "Files of Type" popup.

    (b) In Project/Project Settings, add the "vl\include" directory to
    the C/C++ pane under Category: Preprocessor, "Additional include
    directories". You may have to specify a full or relative path,
    depending on where your project is.

    (c) #include either VLf.h or VLd.h or VLfd.h in your source file. To
    use with the debugging library, add the VL_CHECKING symbol to the
    project settings as in (d) above. (If you don't do this, inline
    functions won't have proper bounds checking.)

Documentation

    Double click on the "vl" file in the doc subdirectory for documentation.

Changes

    1.3.2   Integrated factor routines with main library better.
    1.3.1   Fixed bad cut&paste memory problem in SetSize
    1.3     Fixed problems with stack temporaries being passed as non-const
                references by making *= etc. in-class for Vec[fd].  This is a
                temporary workaround. It breaks some sub vector stuff, such as
                sub(v, 2, 4) *= 3.0; The latest version of gcc shipped with
                redhat 7 no longer takes the -fpermissive flag, however, so it's
                necessary. A redesign to get around stack temporaries no longer
                playing well with constness (and thus necessitating separate
                const and non-const helper classes, and making it harder than
                ever to get C++'s one-conversion-only type system to do what it
                needs to do) is in the works -- it will be VL version 2.0.
            added SubVec -=/+=/etc.
            eliminated v *= m;
            removed all asserts checking that the return value of 'new' is
                not zero. The C++ 'standard' has been changed so that new
                instead throws an exception in this situation.
            added option to use memcpy for general vector copies: VL_USE_MEMCPY.
                The situation seems to have changed so that, on Intel processors
                at least, using memcpy is faster than loops for most cases. See
                note in Vec.cc. (ADVANCED: If you want a container type that
                needs to have its constructor called, this should be disabled.)

    1.2     Solve/factor stuff removed from VLf
                who would ever want heavy duty numerical routines in float?
            changed headers to assume vl/VLd.h, vl/Vec2.h etc
            set VLConfig.h on install, so programs linking against it get right
                flags and don't have to do -D__GCC__ etc. (duh.)
            added Normalise() method to vectors... normalise() is now deprecated.
                (The lower-case vector ops should only be math functions.)
            added zero-length assert to norm and Normalise

    1.1.4   minor VLMath.h fixes
            fixed buggy sparse matrix trans()
            added ConvertVec and ConvertMat for converting between float/double.
                (VLMath.h)
            changed headers to have vl/ prefix (unix version only)
            renamed vl_finite to vl_is_finite
            added vl_inf to VLMath.h
            added vl_rand to VLMath.h
            added clamp[ed] for VecN
            added sign()
            added Factor.cc SVD/QR factorization routines
            updated docs to reflect this
            added xform for 2d & 3d, and mixed versions
            changed makefile structure to current system

    1.1.3   Added linear solver documentation, fixed some windows compile 
            issues under VC++, minor fixes to solver code.

    1.1.2   Added sub, first and last functions for SubVec, and (re?)added
            row function. Added submat functions for SubMat. Added comments
            about the TReal, TVec, etc. types to VL.h

    1.1.1   Fixed problem with += operators and gcc under linux.

    1.1     Added oprod operator, irix n32 CC support, changed makefiles & afs
            support, eliminated most g++ warnings.
