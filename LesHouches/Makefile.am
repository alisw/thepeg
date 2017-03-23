mySOURCES = LesHouchesReader.cc LesHouchesFileReader.cc  \
          LesHouchesEventHandler.cc

DOCFILES = LesHouchesReader.h LesHouchesFileReader.h  \
           LesHouchesEventHandler.h LesHouches.h

INCLUDEFILES = $(DOCFILES) LesHouchesReader.fh \
               LesHouchesFileReader.fh \
                LesHouchesEventHandler.fh

pkglib_LTLIBRARIES = LesHouches.la MadGraphReader.la

# Version info should be updated if any interface or persistent I/O
# function is changed
MadGraphReader_la_LDFLAGS = $(AM_LDFLAGS) -module $(LIBTOOLVERSIONINFO)
MadGraphReader_la_SOURCES = \
MadGraphReader.cc MadGraphReader.h \
MadGraphOneCut.cc MadGraphOneCut.h \
MadGraphTwoCut.cc MadGraphTwoCut.h

# Version info should be updated if any interface or persistent I/O
# function is changed
LesHouches_la_LDFLAGS = $(AM_LDFLAGS) -module -version-info 22:0:0
LesHouches_la_SOURCES = $(mySOURCES) $(INCLUDEFILES)

include $(top_srcdir)/Config/Makefile.aminclude

