// _____________________________________________________________
/* This app is used to plot the output of an e4nu analysis    */
// _____________________________________________________________

#include "plotting/AcceptanceUtils.h"
#include "plotting/XSecUtils.h"
#include <iomanip>
#include <filesystem>

using namespace std; 
using namespace e4nu;
using namespace e4nu::plotting;

/////////////////////////////////////////////////////////////////
// Options:                                                    //
/////////////////////////////////////////////////////////////////
// 1) MC files location                                        //
// 2) Data files location (*)                                  //
// 3) OutputLocation                                           //
// 4) Input MC Files                                           //
// 5) Input data file (*)                                      // 
// 6) Observables List                                         //
// 7) Input Model Names (*)                                    //
// 8) Title (*)                                                //
// 9) Data Name (*)                                            //
// 10) Compute systematics (* def false )                      //
// 11) AddSystematics                                          //
// 12) NoFSI ROOT file                                         //
/////////////////////////////////////////////////////////////////

string mc_location="", data_location="", output_location ="", output_name ="";
vector<string> mc_files, observables, model_names ;
string data_file ="", nofsi_file="", title="", data_name="" ;
bool compute_systematics ;
map<string,double> systematic_map ; 
bool plot_data = true ; 

void PrintFormat(string s);
string GetArg(string op, int argc, char ** argv );
bool ExistArg(string op, int argc, char ** argv );
vector<string> SplitString(string s,char d=',') ; 
int main( int argc, char* argv[] ) {

  std::cout << "Plotting e4nu analysis..." << std::endl;

  if( argc <= 1 ) { 
    PrintFormat("");
    return 0 ; 
  } else { 
    //Loading configuration
   
    if( ExistArg("mc_location",argc,argv)) {
      mc_location = GetArg("mc_location",argc,argv) ;
      cout << "Reading MC files from " << mc_location << std::endl; 
    } else PrintFormat("mc_location") ; 
      
    if( ExistArg("data_location",argc,argv)) {
      data_location = GetArg("data_location",argc,argv) ;
    } else plot_data = false ; 
      
    if( ExistArg("output_location",argc,argv)) {
      output_location = GetArg("output_location",argc,argv) ;

      if( ! std::filesystem::exists(output_location) ) std::filesystem::create_directory(output_location);
    } else PrintFormat("output_location") ; 

    if( ExistArg("output_name",argc,argv)) {
      output_name = GetArg("output_name",argc,argv) ;
    } else PrintFormat("output_name") ; 
      
    string files ; 
    if( ExistArg("input_mc_files",argc,argv)) {
      files = GetArg("input_mc_files",argc,argv) ;
      mc_files = SplitString(files);
    } else PrintFormat("input_mc_files") ; 

    if( ExistArg("input_data_file",argc,argv)) {
      data_file = GetArg("input_data_file",argc,argv) ;
      cout << "Reading data file from " << data_location << data_file << std::endl; 
    } else plot_data = false ; 

    string obs ; 
    if( ExistArg("observable_list",argc,argv)) {
      obs = GetArg("observable_list",argc,argv) ;
      observables = SplitString(obs);
      cout << "Plotting xsec for the following observables: \n- ";
      for( unsigned s = 0 ; s<observables.size();++s ) {
	cout << observables[s] << std::endl;
      }
    } else PrintFormat("observable_list") ; 

    string mdl ; 
    if( ExistArg("model_names",argc,argv)) {
      mdl = GetArg("model_names",argc,argv) ;
      model_names = SplitString(mdl);
      if( model_names.size() != mc_files.size() ){
	std::cout << "Number of mc files does not match the number of models"<< std::endl;
	return 0;
      }
    } 

    cout<<"Loading MC Files:"<<endl;
    for( unsigned int i = 0 ; i < mc_files.size() ; ++i ) { 
      cout << " -> " ;
      if( model_names.size() != 0 ) cout << model_names[i] << ": " ;
      cout << mc_files[i] << endl;
    }
    
    if( ExistArg("title",argc,argv)) {
      title = GetArg("title",argc,argv) ;
    }

    if( ExistArg("data_name",argc,argv)) {
      data_name = GetArg("data_name",argc,argv) ;
    }

    if( ExistArg("nofsi_file",argc,argv)) {
      nofsi_file = GetArg("nofsi_file",argc,argv) ;
      cout << "Adding No FSI file: " << nofsi_file << endl;
    }

    string sys ; 
    if( ExistArg("add-systematics",argc,argv)) {
      sys = GetArg("add-systematics",argc,argv) ;
      vector<string> sys_names = SplitString(sys,':');
      for( unsigned s = 0 ; s < sys_names.size() ; ++s ) { 
	std::cout << sys_names[s]<<std::endl;
	vector<string> tmpsys = SplitString(sys_names[s]) ; 
	if( tmpsys.size() != 2 ) continue ; 
	systematic_map[tmpsys[0]] = stod(tmpsys[1]);
      }
    } 
	
  }

  // Loop over observables
  for( unsigned int i = 0 ; i < observables.size(); ++i ){
    string acceptance_file = plotting::ComputeAcceptance( mc_files, observables[i], title, mc_location, output_location, output_name ) ;
    if( nofsi_file != "" && i == 0 ) { mc_files.push_back(nofsi_file); model_names.push_back("No FSI");}
    plotting::Plot1DXSec( mc_files, data_file, acceptance_file, observables[i], title, data_name, model_names, mc_location, data_location, output_location, output_name, plot_data ) ;
  }

  return 0 ; 
}

void PrintFormat(string s){
  if( s!="") cout << " Missing " << s << endl;
  cout << "RequiredArguments:\n";
  cout << " plote4nuanalysis --mc_location <mc_location> \n --data_location <data_location> \n --output_location <output_loc> \n --output_name <outputname> \n ";
  cout << "--input_mc_files <file1,file2,...,fileN> \n --input_data_file <data> --observable_list <obs1,obs2,...,obsM> " << endl;
  cout << " optional arguments are : \n --model_names <name1,name2,...,nameN> \n --title <title> \n --data_name <data> \n --systematics \n";
  cout << " --add-systematics name,value:name2,value2:...:nameK,valueK \n --nofsi_file <rootfile>" << endl;
}

vector<string> SplitString(string s, char d ) {
  vector<string> strings;
  int startIndex = 0, endIndex = 0;
  for (unsigned int i = 0; i <= s.size(); i++) {
        
    // If we reached the end of the word or the end of the input.
    if (s[i] ==d || i == s.size()) {
      endIndex = i;
      string temp;
      temp.append(s, startIndex, endIndex - startIndex);
      strings.push_back(temp);
      startIndex = endIndex + 1;
    }
  }
  return strings ;   
}

string GetArg(string op, int argc, char ** argv )
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  while(argc>2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-') {

	char op_cur[buf_size];
	strcpy(op_cur,&argv[1][2]);

	if (strcmp(op.c_str(),op_cur)==0) {
	  if (strlen(&argv[2][0]) ) {
            strcpy(argument,&argv[2][0]);
	  } 
	}
      }
      argc--;
      argv++;

    }

  string value = string(argument);
  delete [] argument;
  return value ; 
}


bool ExistArg(string op, int argc, char ** argv )
{
  const int buf_size = 2048*128;
  char *  argument   = new char[buf_size];
  strcpy(argument, "");

  while(argc>2)
    {
      if (argv[1][0] == '-' && argv[1][1] == '-') {

	char op_cur[buf_size];
	strcpy(op_cur,&argv[1][2]);

	if (strcmp(op.c_str(),op_cur)==0) {
	  return true ; 
	} 
      }
      argc--;
      argv++;

    }
  delete [] argument ; 
  return false;
}
