CC = g++

FLAG = -O3 -g
RUNFILE = GEB_HFC

OBJFILES = GEB_HFC.o HFC.o

.PHONY: all clean

all: $(RUNFILE)

$(RUNFILE): $(OBJFILES)
	$(CC) $(FLAG) $(OBJFILES) -o $(RUNFILE)

GEB_HFC.o: GEB_HFC.cpp HFC.h ntoh64.h global.h
	$(CC) $(FLAG) -c GEB_HFC.cpp

HFC.o: HFC.cpp HFC.h ntoh64.h global.h
	$(CC) $(FLAG) -c HFC.cpp

clean:
	rm -f $(RUNFILE) $(OBJFILES)
