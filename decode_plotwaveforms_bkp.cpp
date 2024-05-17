#include <stdio.h>
#include <iostream>
#include <iterator>
#include "TFile.h"
#include <TTree.h>
#include <TCanvas.h>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <TF1.h>
#include <vector>
#include <cstdlib>
#include <filesystem>
using namespace std;



// a long long int can store 8 bytes (64 bits) of data
typedef int BLKSIZE;

// get file size in bytes
long getFileSize(FILE * file) {
  long lCurPos, lEndPos;
  lCurPos = ftell(file);
  fseek(file, 0, 2);
  lEndPos = ftell(file);
  fseek(file, lCurPos, 0);
  return lEndPos;
}

// Function to remove ".csv" extension from the file path
std::string removeCsvExtension(const std::string& filePath) {
    if (filePath.size() > 4 && filePath.substr(filePath.size() - 4) == ".csv") {
        return filePath.substr(0, filePath.size() - 4);
    }
    return filePath;
}

// main routine
int decode_plotwaveforms_st(std::string outdir, std::vector<int> file_numbers) {

  const int datacheck = 0; //enter 1 if checking fake data, 2 if checking baseline, 0 for no checking
  //will check every 20th adc value
  const int dataplot = 1; //enter 1 if you want the waveforms plotted
  // vector<int> file_numbers;
  // cout << "Enter the file numbers you'd like to decode (separated by spaces): ";
  // string line;
  // getline(cin, line);
  // istringstream iss(line);
  // file_numbers = vector<int>(istream_iterator<int>(iss), istream_iterator<int>());

  double maxamps[file_numbers.size()];
  double ey[file_numbers.size()];

  //path definitions
  const char *path0 = "./data/";
  const char *path1 = "xmit_trig_bin_grams_";
  const char *path3 = ".dat";
  const char *opath1 = "outfile_";
  const char *opath2 = ".root";

  for(vector<int>::iterator it = file_numbers.begin(); it != file_numbers.end(); ++it){
    // convert to string
    int f = *it;
    std::stringstream ss;
    ss << f;
    std::string file_number = ss.str();


    std::string fullpath = std::string(path0) + std::string(path1) + file_number + std::string(path3);
    std::string ofullpath = std::string(outdir) + std::string(opath1) + file_number + std::string(opath2);
    
    const char *filePath = fullpath.c_str();
    BLKSIZE *fileBuf; // pointer to buffered data
    FILE *file = NULL; // file pointer 
    long fileSize;

    int debug = 0; //set to 1 for more verbose output

    // some decoding variables
    int hdcount = 0; //keeps track of how many header words found in channel
    int chcount = 64; //keeps track of how many channels found in FEM
    int femcount = 0; //keeps track of how many fems found in event
    int evtcount = 0; //keeps track of how many events found in all files read
    int thisevtcount = 0; //keep tracks of how many events found in each file
    int adccount = 0; //keeps track of how many adc words were read in
    int samplecount = 0; //keeps track of how many adc values were read in
    int badword;
    // variables which are part of FEM header
    int hmodadd;
    int hmodid;
    int hnwords; //number of adc words to read
    int hevt;
    int hframe;
    int hchecksum;
    int htrigsamp;
    int htrigfram;

    //variables which are part of channel header
    int ch=0;

    static const int FEMs=20;
    static const int CHs =64;
    static const int MAXsamp=250; 
    const char * ofilePath = ofullpath.c_str();

    //channel vs ADC val histograms
    TFile* out_file = new TFile(ofilePath,"RECREATE");

    //create TTree
    TTree* adc_tree = new TTree("adc_tree", "ADC Counts");
    int fem_id, channel_id, event_id, adc_count;
    adc_tree->Branch("fem_id", &fem_id, "fem_id/I");
    adc_tree->Branch("channel_id", &channel_id, "channel_id/I");
    adc_tree->Branch("event_id", &event_id, "event_id/I");
    adc_tree->Branch("adc_count", &adc_count, "adc_count/I");


    int i, e, d;

    //process binary files one by one
    //  for (int e=0; e<1; e++){ //run over files
    //    for (int d=0; d<1; d++){ //run over files
    //      sprintf(filePath,"./data/xmit_trig_bin_grams_100.dat");//%i_%i.dat",e,d);//input file

    // open binary file and check if open was successful

  	printf("analysing file%s",filePath);

    if ((file = fopen(filePath, "rb")) == NULL) {
      cout << "\nCould not open file:\t" << filePath << endl;
    } else {
      if (debug == 1) cout << "\nFile\t" << filePath << "\topened successfully." << endl;

      // get file size in bytes
      fileSize = getFileSize(file);
      if (debug == 1) cout << "File size is " << fileSize << " bytes." << endl;
      if (debug == 1) cout << "Block size is " << sizeof(BLKSIZE) << " bytes." << endl;

      // allocate space in the buffer for the whole file
      fileBuf = new BLKSIZE[fileSize / sizeof(BLKSIZE)];

      // read whole file into buffer
      fread(fileBuf, fileSize, 1, file);

      thisevtcount = 0;

      //=========================================================
      // decode, one BLKSIZE-length buffer at a time
      for (i = 0; i < fileSize / sizeof(BLKSIZE); i++) {

        //decode 32-bit word
        badword=0;
        if ((fileBuf[i] & 0xffffffff) == 0xffffffff) { //begin of event
          printf("\n=====> New event %i\n", evtcount);
          femcount = 0;
          evtcount++;
          thisevtcount++;
        } else if ((fileBuf[i] & 0xffffffff) == 0xe0000000) { //end of event
          printf("=====> End of event\n");
          if (femcount != 2) {
            printf("############# ERROR: Missing FEMs! #############\n");
            printf("              Found %d FEMs instead of 2\n", femcount);
          }
        } else { //split in two 16-bit words and decode

          //check lower 16 bits first

          if (((fileBuf[i] & 0xffff) & 0xf000) == 0xf000) { //if header word
            if (debug == 1) printf("Debug: -----> FEM header word\n");
            if ((((fileBuf[i] & 0xffff) & 0xffff) == 0xffff) && (hdcount==12 || hdcount==0)) { //first header word
              if (chcount != 64) {
                printf("############# ERROR: Missing channels! #############\n");
                printf("              Found %d channels instead of 64\n", chcount);
              }
              if (adccount != hnwords && (adccount - 1) != hnwords) {
                printf("############# ERROR: Wrong number of adc words! #############\n");
                printf("                     Read %i, Expected from header %i(+1)\n", adccount, hnwords);
              }
              hdcount = 1;
              femcount++;
              chcount = 0;
              adccount = 0;
            } //end if first header word
            else { //following header words
              hdcount++;
  	    printf("hdcount : %i\n",hdcount);
              if (hdcount == 2) {
                hmodadd = ((fileBuf[i] >> 16) & 0xfff) & 0x1f;
                hmodid = (((fileBuf[i] >> 16) & 0xfff) >> 5) & 0x7f;
                printf("\t---\n");
                printf("\tModule address: %d\n", hmodadd);
              } else if (hdcount == 3) {
                hnwords = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
              } else if (hdcount == 4) {
                //part of hnwords, above
              } else if (hdcount == 5) {
                //if (hevt == 5) break;
                hevt = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
                printf("\tEvent number  : %d\n", hevt);

              } else if (hdcount == 6) {
                //part of hevt, above
              } else if (hdcount == 7) {
                hframe = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
                printf("\tFrame number  : %d\n", hframe);
              } else if (hdcount == 8) {
                //part of frame number, above
              } else if (hdcount == 9) {
                hchecksum = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
              } else if (hdcount == 10) {
                //part of checksum, above
              } else if (hdcount == 11) {
                htrigfram = (fileBuf[i]) & 0xff; //lower 4 bits of trigger frame number [7:4] and upper 4 bits of 2MHz trigger sample number [3:0] (dictates when trig was received by FEM)
                htrigsamp = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xff) << 12); //lower 8 bits of trigger 2MHz sample number, plus previous 4 bits (dictates when trig was received by FEM)
  	      printf("htrigfram : %d, htrigsamp : %d\n",htrigfram, htrigsamp);
              } else if (hdcount == 12) {
                //part of trigtime, above
              } else {
                printf("############# ERROR: Unknown header word! #############\n");
                printf("              %x\n", (fileBuf[i]) & 0xffff);
              }
            } //end if following header
          } //end if header word
          else if (((fileBuf[i] & 0xffff) & 0xf000) == 0x4000) { //channel begin
            if (hdcount != 12) {
              printf("############# ERROR: Missing headers! #############\n");
              printf("              Found %d headers instead of 12\n", hdcount);
            }
            if (debug == 1) printf("Debug: Channel begin\t");
            ch = fileBuf[i] & 0xfff;
            chcount++;
            adccount++;
            samplecount = 0;
          } else if (((fileBuf[i] & 0xffff) & 0xf000) == 0x5000) { //channel end
            if (debug == 1) printf("Debug: Channel end\t(%d)\n", ch);
            adccount++;
          } else {
            adccount++;
            samplecount++;
            //fill adc values for channel 'ch', module 'hmodadd', and event 'hevt' here
            //if (dataplot == 1) waveform[hmodadd][ch][hevt] -> SetBinContent(samplecount, fileBuf[i] & 0xfff); //12 bits
            fem_id = hmodadd;
            channel_id = ch;
            event_id = hevt;
            adc_count = (fileBuf[i] & 0xfff);
            adc_tree->Fill();
            //printf("ADC value %d \n",fileBuf[i]&0xfff);
            if (datacheck == 1 && (samplecount % 20 == 0)) { //check fake data is correct 
              if ((fileBuf[i] & 0xffff) != 0xa33) { //a33 is the fake data pattern 
                printf("############# ERROR: Found bad ADC data! #############\n");
                cout << "\tFile: " << filePath << endl;
                printf("Event\tFEM\tChannel\tSample\tADC value\n");
                printf("%d\t%d\t%d\t%d\t%x\n", hevt, hmodadd, ch, samplecount, (fileBuf[i] & 0xffff));
              }
            } else if (datacheck == 2 && (samplecount % 20 == 0)) {
              if ((fileBuf[i] & 0xff00) != 0x800 && (fileBuf[i] & 0xff00) != 0x700 && (fileBuf[i] & 0xff00) != 0x100) { //baseline 
                printf("############# ERROR: Found bad ADC data! #############\n");
                cout << "\tFile: " << filePath << endl;
                printf("Event\tFEM\tChannel\tSample\tADC value\n");
                printf("%d\t%d\t%d\t%d\t%x\n", hevt, hmodadd, ch, samplecount, fileBuf[i] & 0xffff);
              }
            }
          } //end if adc value word 

          //check upper 16 bits next
          if ((((fileBuf[i] >> 16) & 0xffff) & 0xf000) == 0xf000) { //header word
            if (debug == 1) printf("Debug: -----> FEM header word\n");
            if (((((fileBuf[i] >> 16) & 0xffff) & 0xffff) == 0xffff) && (hdcount==12 || hdcount==0))  { //first header word
              if (chcount != 64) {
                printf("############# ERROR: Missing channels! #############\n");
                printf("              Found %d channels instead of 64\n", chcount);
              }
              if (adccount != hnwords && (adccount - 1) != hnwords) {
                printf("############# ERROR: Wrong number of adc words! #############\n");
                printf("                     Read %i, Expected from header %i(+1)\n", adccount, hnwords);
              }
              hdcount = 1;
              femcount++;
              chcount = 0;
              adccount = 0;
            } //end if first header word
            else { //following header words
              hdcount++;
              if (hdcount == 2) {
                hmodadd = ((fileBuf[i] >> 16) & 0xfff) & 0x1f;
                hmodid = (((fileBuf[i] >> 16) & 0xfff) >> 5) & 0x7f;
                printf("\t---\n");
                printf("\tModule address: %d\n", hmodadd);
              } else if (hdcount == 3) {
                hnwords = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
              } else if (hdcount == 4) {
                //part of hnwords, above
              } else if (hdcount == 5) {
                hevt = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
                printf("\tEvent number  : %d\n", hevt);
              } else if (hdcount == 6) {
                //part of hevt, above
              } else if (hdcount == 7) {
                hframe = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
                printf("\tFrame number  : %d\n", hframe);
              } else if (hdcount == 8) {
                //part of frame number, above
              } else if (hdcount == 9) {
                hchecksum = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xfff) << 12);
              } else if (hdcount == 10) {
                //part of checksum, above
              } else if (hdcount == 11) {
                htrigfram = (fileBuf[i]) & 0xff; //lower 4 bits of trigger frame number [7:4] and upper 4 bits of 2MHz trigger sample number [3:0] (dictates when trig was received by FEM)
                htrigsamp = ((fileBuf[i] >> 16) & 0xfff) + ((fileBuf[i] & 0xff) << 12); //lower 8 bits of trigger 2MHz sample number, plus previous 4 bits (dictates when trig was received by FEM)
              } else if (hdcount == 12) {
                //part of trigtime, above
              } else {
                printf("############# ERROR: Unknown header word! #############\n");
                printf("              %x\n", (fileBuf[i] >> 16) & 0xffff);
                //badword=1;
              }
            } //end if following header words
          } //end if header word
          else if ((((fileBuf[i] >> 16) & 0xffff) & 0xf000) == 0x4000) { //channel begin
            if (hdcount != 12) {
              printf("############# ERROR: Missing headers! #############\n");
              printf("              Found %d headers instead of 12\n", hdcount);
              //badword=1;
            }
            if (debug == 1) printf("Debug: Channel begin\t");
            ch = (fileBuf[i] >> 16) & 0xfff;
            chcount++;
            adccount++;
            samplecount = 0;
          } else if ((((fileBuf[i] >> 16) & 0xffff) & 0xf000) == 0x5000) { //channel end
            if (debug == 1) printf("Debug: Channel end\t(%d)\n", ch);
            adccount++;
          } else { //adc value 
            adccount++;
            samplecount++;
            //fill adc values for channel 'ch', module 'hmodadd', and event 'hevt' here
            //if (dataplot == 1) waveform[hmodadd][ch][hevt] -> SetBinContent(samplecount, (fileBuf[i] >> 16) & 0xfff); //12 bits
            fem_id = hmodadd;
            channel_id = ch;
            event_id = hevt;
            //if(badword=0){
            adc_count = (fileBuf[i] >> 16) & 0xfff;
            adc_tree->Fill();
            //}
            //printf("ADC value %d \n",fileBuf[i]&0xfff);
            if (datacheck == 1 && (samplecount % 20 == 0)) { //check fake data is correct
              if (((fileBuf[i] >> 16) & 0xffff) != 0xa33) { //a33 is the fake data pattern 
                printf("############# ERROR: Found bad ADC data! #############\n");
                cout << "\tFile: " << filePath << endl;
                printf("Event\tFEM\tChannel\tSample\tADC value\n");
                printf("%d\t%d\t%d\t%d\t%x\n", hevt, hmodadd, ch, samplecount, (fileBuf[i] >> 16) & 0xffff);
              }
            } else if (datacheck == 2 && (samplecount % 20 == 0)) {
              if (((fileBuf[i] >> 16) & 0xff00) != 0x800 && ((fileBuf[i] >> 16) & 0xff00) != 0x700 && ((fileBuf[i] >> 16) & 0xff00) != 0x100) { //baseline 
                printf("############# ERROR: Found bad ADC data! #############\n");
                cout << "\tFile: " << filePath << endl;
                printf("Event\tFEM\tChannel\tSample\tADC value\n");
                printf("%d\t%d\t%d\t%d\t%x\n", hevt, hmodadd, ch, samplecount, (fileBuf[i] >> 16) & 0xffff);
              }
            }
          } //end if adc value word 

        } //end else (split in two words and decode)

      } //end for loop

      if (debug == 1) printf("Found %i new events so far; %i new events in this buffer.\n", evtcount, thisevtcount);

      if ((fileBuf[fileSize / sizeof(BLKSIZE) - 1] & 0xffffffff) != 0xe0000000) {
        if (debug == 1) printf("%%WARNING: buffer does not end on 'last event word'. Event may be incomplete.\n");
      }

      delete fileBuf;
      fclose(file);
      if (debug == 1) cout << "File closed. Buffer deleted.\n" << endl;

    } //end if file exists
    //    }//end loop over dma ops
    //}//end loop over events to read
    adc_tree->Write();
    out_file->Close();	
  }
return 0;

}

int main(int argc, char* argv[]) {

  std::string csvFilePath = argv[1];
  std::ifstream csvFile(csvFilePath);
  std::string outdir = removeCsvExtension(csvFilePath);

  if (!csvFile.is_open()) {
      std::cerr << "Error: Could not open file " << csvFilePath << std::endl;
      return 1;
      }

  std::string line;
  std::vector<int> file_numbers;

  // Skip the header line
  if (std::getline(csvFile, line)) {
  }

  // Read the file line by line
  while (std::getline(csvFile, line)) {
      std::stringstream ss(line);
      std::string cell;

      if (std::getline(ss, cell, ',')) {
        int fileId = std::stoi(cell);
        file_numbers.push_back(fileId);
      }
  }

  csvFile.close();

  //mkdir(outdir);

  if (!std::filesystem::exists(outdir)) {
      if (!std::filesystem::create_directories(outdir)) {
          std::cerr << "Error: Could not create directory " << outdir << std::endl;
          return 1;
      }
  }

  decode_plotwaveforms_st(outdir, file_numbers);

  return 0;

}
