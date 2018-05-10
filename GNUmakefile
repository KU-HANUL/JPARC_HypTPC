name := main
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

#root
CPPFLAGS += `root-config --cflags`
LDFLAGS  += `root-config --libs`

#event gen
EVTGENDIR       = /home/sbyang/software/EvtGen/EvtGen/EvtGen/R01-01-00
EVTGENINCDIR    = $(EVTGENDIR)
EVTGENLIBDIR    = $(EVTGENDIR)/lib
EVTGENINCLUDE   = -I$(EVTGENINCDIR)
EVTGENLIB       = -L$(EVTGENLIBDIR) -lEvtGen -lEvtGenExternal
EVTGEN_EXTERNAL = 0
EXTRAFLAGS      = -DEVTGEN_EXTERNAL

CPPFLAGS += $(EVTGENINCLUDE)
LDFLAGS += $(EVTGENLIB)

#HEPMC
HEPMCDIR       = /home/sbyang/software/EvtGen/HepMC/2.06.09
HEPMCINCDIR    = $(HEPMCDIR)/include
HEPMCLIBDIR    = $(HEPMCDIR)/lib
HEPMCINCLUDE   = -I$(HEPMCINCDIR)
HEPMCLIB       = -L$(HEPMCLIBDIR) -lHepMC

CPPFLAGS += $(HEPMCINCLUDE)
LDFLAGS += $(HEPMCLIB)

#Pythia
PYTHIADIR       = /home/sbyang/software/EvtGen/PYTHIA/pythia8170
PYTHIAINCDIR    = $(PYTHIADIR)/include
PYTHIALIBDIR    = $(PYTHIADIR)/lib
PYTHIAINCLUDE   = -I$(PYTHIAINCDIR)
PYTHIALIB       = -L$(PYTHIALIBDIR) -lpythia8 -llhapdfdummy -lhepmcinterface

CPPFLAGS += $(PYTHIAINCLUDE)
LDFLAGS += $(PYTHIALIB)

#PHOTOS                                                                                                     
PHOTOSDIR       = /home/sbyang/software/EvtGen/PHOTOS/PHOTOS
PHOTOSINCDIR    = $(PHOTOSDIR)/include
PHOTOSLIBDIR    = $(PHOTOSDIR)/lib
PHOTOSINCLUDE   = -I$(PHOTOSINCDIR)
PHOTOSLIB       = -L$(PHOTOSLIBDIR) -lPhotosCxxInterface -lPhotosFortran

CPPFLAGS += $(PHOTOSINCLUDE)
LDFLAGS += $(PHOTOSLIB)


#CLHEP
CLHEPDIR        = /home/sbyang/software/CLHEP/clhep-2.3.3.0-install
CLHEPINCDIR     = $(CLHEPDIR)/include/CLHEP
CLHEPLIBDIR     = $(CLHEPDIR)/lib
CLHEPINCLUDE    = -I$(CLHEPINCDIR)
CLHEPLIB        = -L$(CLHEPLIBDIR) -lCLHEP

CPPFLAGS += $(CLHEPINCLUDE)
LDFLAGS += $(CLHEPLIB)

#X11
#X11DIR         = /usr/X11R6                                                                      
X11DIR          = /usr
X11LIBDIR       = $(X11DIR)/lib
X11INCDIR       = $(X11DIR)/include
XINCLUDE        = -I$(X11INCDIR)
XLIBS           = -L$(X11LIBDIR) -lX11 -lXmu

CPPFLAGS += $(XINCLUDE)
LDFLAGS += $(XLIB)
