# RFC 2: CostalME Committer Guidelines

Author: Andres Payo

Contact: agarcia at bgs.ac.uk

Status: Draft

Purpose
-------

CoastalME works with Git, a distributed version control system offering full-scale developer support via GitHub cloud-based hosting system. To formalize Source Version Control (or SVC) commit access, and specify some guidelines for SVC committers.

Election to SVC Commit Access
-----------------------------

Permission for SVC commit access shall be provided to new developers
only if accepted by the CoastalME Project Steering Committee. A proposal
should be written to the PSC for new committers and voted on normally.
It is not necessary to write an RFC document for these votes, a proposal
to [CoastalME Discussions forum](https://github.com/apayo/CoastalME/discussions/19#discussion-7445785) is sufficient.

Removal of SVC commit access should be handled by the same process.

The new committer should have demonstrated commitment to CoastalME and
knowledge of the CoastalME source code and processes to the committee's
satisfaction, usually by reporting bugs, submitting patches, and/or
actively participating in the CoastalME mailing list(s).

The new committer should also be prepared to support any new feature or
changes that he/she commits to the CoastalME source tree in future
releases, or to find someone to which to delegate responsibility for
them if he/she stops being available to support the portions of code
that he/she is responsible for.

All committers should also be a member of the coastalme-dev mailing list so
they can stay informed on policies, technical developments and release
preparation.

New committers are responsible for having read, and understood this
document.

Committer Tracking
------------------

A list of all project committers will be kept in the main CoastalME GitHub directory
(called [COMMITTERS](../../COMMITERS.md)) listing for each SVC committer:

-  Userid: the id that will appear in the SVC logs for this person.
-  Full name: the users actual name.
-  Email address: A current email address at which the committer can be
   reached. It may be altered in normal ways to make it harder to
   auto-harvest.
-  A brief indication of areas of responsibility.

SVC Administrator
-----------------

One member of the Project Steering Committee will be designed the SVC
Administrator. That person will be responsible for giving SVC commit
access to folks, updating the COMMITTERS file, and other SVC related
management. That person will need login access on the SVC server of
course.

Initially David Favis-Mortlock will be the SVC Administrator.

SVC Commit Practices
--------------------

The following are considered good SVC commit practices for the CoastalME
project.

-  Use meaningful descriptions for SVC commit log entries.
-  Add a bug reference like "(#1232)" at the end of SVC commit log
   entries when committing changes related to a ticket in Trac. The '#'
   character enables Trac to create a hyperlink from the changeset to
   the mentioned ticket.
-  After committing changes related to a ticket in Trac, write the tree
   and revision in which it was fixed in the ticket description. Such as
   "Fixed in trunk (r12345) and in branches/1.7 (r12346)". The 'r'
   character enables Trac to create a hyperlink from the ticket to the
   changeset.
-  Changes should not be committed in stable branches without a
   corresponding bug id. Any change worth pushing into the stable
   version is worth a bug entry.
-  Never commit new features to a stable branch without permission of
   the PSC or release manager. Normally only fixes should go into stable
   branches.
-  New features go in the main development trunk.
-  Only bug fixes should be committed to the code during pre-release
   code freeze, without permission from the PSC or release manager.
-  Significant changes to the main development version should be
   discussed on the gdal-dev list before you make them, and larger
   changes will require a RFC approved by the PSC.
-  Do not create new branches without the approval of the PSC. Release
   managers are assumed to have permission to create a branch.
-  All source code in SVC should be in Unix text format as opposed to
   DOS text mode.
-  Please make use of existing code (such as CoastalME's 2d integer and floating point classes) wherever possible. Don't try to re-invent the wheel!
-  When committing new features or significant changes to existing
   source code, the committer should take reasonable measures to ensure
   that the source code continues to build and work on the most commonly
   supported platforms (currently Linux, Windows, and Mac), either by testing
   on those platforms directly, running [wiki:Buildbot] tests, or by
   getting help from other developers working on those platforms. If new
   files or library dependencies are added, then the configure.in,
   Makefile.in, Makefile.vc and related documentations should be kept up
   to date.
   
CoastalME coding guidelines
---------------------------

- Clarity, readability, and consitency of style are the most important attributes of any source code submitted to the CoastalME GitHub repository. The easy understandability of newly-submitted code is much more important than its run-time efficiency. Code can always be optimized later. "Premature optimization is the root of all evil", according to Donald Knuth (https://en.wikipedia.org/wiki/Donald_Knuth) 

- If you submit source code which isn't clear, readable, well-commented, and in a style similar to existing CoastalME source code, then we'll ask you to tidy it then re-submit. We hope that you don't mind.

- Here are some guidelines:

   - Do not save CoastalME source code files (*.cpp, *.h) with a limited line length (e.g. 80 characters). One consequence of doing this is that headers for fixed-width output (e.g. PER_ITER_HEAD1) are split and redistributed over several lines. This makes it much harder to ensure correct run-time alignment of such headers. Instead, allow unlimited line length. If this is not possible, allow a line length of at least 250 characters.

   - Please name your variables using (as far as possible) Hungarian Notation (https://en.wikipedia.org/wiki/Hungarian_notation)
   
   - Our preferred coding style is very broadly based on that of Microsoft (see https://clang.llvm.org/docs/ClangFormatStyleOptions.html#basedonstyle) but with three spaces as indentation. By far the easiest way to find out more about the CoastalME coding style is to look at the existing source code.
   
   - You may find it helpful to use a source code formatter before submitting your code. The clang-format source code formatter does an OK job (it is still under active development). To install this on a Linux system, run "apt install clang-format" or equivalent. There should be a clang-format configuration file (.clang-format) in your src directory if you have pulled the CoastalME code from GitHub. To use clang-format to tidy CoastalME source code, run tidy_src.sh <filename> in the src directory. This will modify <filename> (which can be a .cpp or a .h file) "in place" i.e. without making a backup. Please DO NOT run clang-format on files which you have not modified. For more about clang-format, see https://clang.llvm.org/docs/ClangFormat.html

   - Fix all compiler warnings relating to your code before submitting it.

   - Also you are strongly advised to use a static code analyzer such as clang-tidy on your new or modified source code. Run do_clang-analyse.sh in the src directory, this will produce a text file named 000_clang-analyze_advice.txt, which gives suggestions for e.g. maiaining const-correctness of variables. Note however that not all suggestions by clang-tidy "work", some may prevent your code building. For this reason, it is sensible to try just a few suggestions at a time, then do a test build. Then repeat, trying more clang-analyse suggestions.
   
   - It is also a good idea to try your modified version of CoastalME with the example datasets which are on GitHub (this is what we will do, to test your contribution, so you may as well try it first!). Use run_test_suite.sh for this. 
   

Relationship with other upstream projects imported in CoastalME code base
------------------------------------------------------------------------

Some parts of the CoastalME code base are regularly refreshed from other
upstream projects. So changes in those areas should go first into those
upstream projects, otherwise they may be lost during a later refresh.

Currently the list of those areas is :

-  https://github.com/OSGeo/gdal/
-  https://github.com/erdc/cshore

Legal
-----

Committers are the front line gatekeepers to keep the code base clear of
improperly contributed code. It is important to the CoastalME users,
developers and the OSGeo foundation to avoid contributing any code to
the project without it being clearly licensed under the project license.

Generally speaking the key issues are that those providing code to be
included in the repository understand that the code will be released
under the GNU license, and that the person providing the code has the
right to contribute the code. For the committer themselves understanding
about the license is hopefully clear. For other contributors, the
committer should verify the understanding unless the committer is very
comfortable that the contributor understands the license (for instance
frequent contributors).

If the contribution was developed on behalf of an employer (on work
time, as part of a work project, etc) then it is important that an
appropriate representative of the employer understand that the code will
be contributed under the GNU license. The arrangement should be
cleared with an authorized supervisor/manager, etc.

The code should be developed by the contributor, or the code should be
from a source which can be rightfully contributed such as from the
public domain, or from an open source project under a compatible
license.

All unusual situations need to be discussed and/or documented.

Committers should adhere to the following guidelines, and may be
personally legally liable for improperly contributing code to the source
repository:

-  Make sure the contributor (and possibly employer) is aware of the
   contribution terms.
-  Code coming from a source other than the contributor (such as adapted
   from another project) should be clearly marked as to the original
   source, copyright holders, license terms and so forth. This
   information can be in the file headers, but should also be added to
   the project licensing file if not exactly matching normal project
   licensing [LICENSE.md](../../LICENSE.md).
-  Existing copyright headers and license text should never be stripped
   from a file. If a copyright holder wishes to give up copyright they
   must do so in writing to the foundation before copyright messages are
   removed. If license terms are changed it has to be by agreement
   (written in email is ok) of the copyright holders.
-  Code with licenses requiring credit, or disclosure to users should be
   added to [LICENSE.md](../../LICENSE.md).
-  When substantial contributions are added to a file (such as
   substantial patches) the author/contributor should be added to the
   list of copyright holders for the file.
-  If there is uncertainty about whether a change it proper to
   contribute to the code base, please seek more information from the
   project steering committee, or the foundation legal counsel.

Bootstrapping
------------

The following existing committers will be considered authorized CoastalME
committers as long as they each review the committer guidelines, and
agree to adhere to them. The SVC administrator will be responsible for
checking with each person.

-  David Favis-Mortlock
-  Andres Payo
-  Manuel Cobos Budia

--------------

-  [COMMITERS](../../COMMITERS.md)
   
