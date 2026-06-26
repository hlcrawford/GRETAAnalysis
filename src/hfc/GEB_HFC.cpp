#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <fcntl.h>

#include <arpa/inet.h>

#include <signal.h>

#include "global.h"
#include "HFC.h"
#include "ntoh64.h"

#define SQR(x)         ((x)*(x))

using namespace std;

int gotsignal;
void breakhandler(int dummy) {
  gotsignal = 1;
}


void swapbytes(char* a, char *b)
{
  char tmp=*a;
  *a=*b;
  *b=tmp;
}


void BrowseData(gebData header) {
  cerr << "type:" << header.type 
       << " len: " << header.length
       << " ts: 0x" << hex << header.timestamp << dec
       << endl;
}

int main(int argc, char** argv) {

  gotsignal = 0;
  signal(SIGINT, breakhandler);
  signal(SIGPIPE, breakhandler);

  if(argc==1) {
    cerr << argv[0] << " <flag: -p (pipeout) or -z (.gz input file)> <Input file>" << endl
	 << "brings GRETINA Mode3 event file" << endl
	 << "in proper sequence" << endl;
    exit(0);
  }
  
  bool pipeflag = false;
  bool zipflag = false;
  bool bzipflag = false;
  FILE *in = NULL;
  FILE *out = NULL;

  string filename;
  while (argc > 1) {
    if (!(strcmp(argv[1], "-p"))) {
      argc--; argv++;
      pipeflag = true;
    } else if (!(strcmp(argv[1], "-z"))) {
      argc--; argv++;
      zipflag = true;
    } else if (!(strcmp(argv[1], "-bz"))) {
      argc--; argv++;
      bzipflag = true;
    } else {
      filename = argv[1];
      argc--; argv++;
    }
  }

  if (pipeflag == true) {
    out = stdout;
  } else {
    out = fopen("HFC.dat", "wb");
  }
  
  if (!zipflag && !bzipflag) {
    in = fopen(filename.c_str(), "rb");
    if (!in) {
      if (!pipeflag) {
	cerr << "HFC: cannot open file " << filename.c_str() << endl;
      }
      return 1;
    } else {
      if (!pipeflag) {
	cout << "HFC: opened file " << filename.c_str() << endl;
      }
    }
  } else if (zipflag && !bzipflag) {
    string zfilename = "zcat " + filename;
    in = popen(zfilename.c_str(), "r");
    if (!in) {
      if (!pipeflag) {
	cerr << "HFC: cannot open file " << zfilename.c_str() << endl;
      }
      return 1;
    } else {
      if (!pipeflag)
	cout << "HFC: opened file " << zfilename.c_str() << endl;
    }
  } else if (bzipflag && !zipflag) {
    string zfilename = "bzcat " + filename;
    in = popen(zfilename.c_str(), "r");
    if (!in) {
      if (!pipeflag) {
	cerr << "HFC: cannot open file " << zfilename.c_str() << endl;
      }
      return 1;
    } else {
      if (!pipeflag)
	cout << "HFC: opened file " << zfilename.c_str() << endl;
    }
  }
  
  long long totread = 0;
  int read;
  int EvtCount=0;
  BYTE cBuf[8*16382];
  gebData aGeb;
  HFC hfc_list(500*8192, out);
  // 972: strange mode 2 with mem depth 40*8192 needed
  
  bool success=true;
	    
  while (fread(&aGeb, sizeof(struct gebData), 1, in)==1 && !gotsignal) {

    aGeb.length = ntohs(aGeb.length);
    aGeb.seqnum = ntohs(aGeb.seqnum);
    aGeb.timestamp = (int64_t)ntoh64((uint64_t)aGeb.timestamp);
    
    read = fread(cBuf, sizeof(char), aGeb.length, in);
    totread += read + sizeof(struct gebData);
    if (read != aGeb.length) {
      if (!pipeflag) {
	cerr << aGeb.length << " bytes expected but"
	     << read << " bytes read. Bailing out"
	     << endl;
	cerr.flush();
      }
      break;
    }
    
    EvtCount++;
		
    if((EvtCount % 10000)==0 && !pipeflag) {
      cerr << "Event " << EvtCount
	   << " read:" << read
	   << " total read:" << totread/1000000
	   << " Mb \r";
      cerr.flush();
    }
    
#if(0)
    cerr << "Event:" << EvtCount 
	 << " #data:" << read
	 << " geb:" << sizeof(struct gebData)
	 << " total bytes read: " << totread
	 << " (0x" << hex << totread << dec << ")"
	 << endl;
#endif	
    
    if(!hfc_list.add(aGeb, cBuf) && success) {
      success = false;
      if (!pipeflag) {
	cerr << "HFC: adding event in HFC failed"
	     << endl;
      }
    }
  }
  
  if (!pipeflag) {
    cerr << "HFC: calling flush" << endl; cerr.flush();
  }
  hfc_list.flush();
  hfc_list.printstatus();
  
  if (!pipeflag) {
    cerr << "HFC: closing files" << endl; cerr.flush();
  }
  if (zipflag || bzipflag)
    pclose(in);
  else
    fclose(in);
  if (!pipeflag)
    fclose(out);
  if (!pipeflag) {
    cerr << "HFC: done" << endl; cerr.flush(); 
  }
}

	 
