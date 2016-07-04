#include <string>
#include "base/CommandLineParser.h"
#include "src/DNAVector.h"
#include "visual/Histogram.h"
#include "base/FileParser.h"
#include "src/DataRead.h"



class Command
{
public:
  Command() {
    Add("basic", "BasicStats");
    Add("report", "MakeBasicReport");
    Add("lapcands", "FindOverlapCands");
    Add("esterr", "EstimateErrors");
    Add("lapstats", "LapStats");

  }

  void Add(const string & cmd, const string & exe, const string & dep = "") {
    m_cmd.push_back(cmd);
    m_program.push_back(exe);
    m_depend.push_back(dep);
  }


  void PrintLogo() const {

    cout << "    _____      _____            _" << endl;
    cout << "   |_   _|    /  ___|          | |" << endl;
    cout << "     | | _ __ \\ `--.  ___  __ _| |_" << endl;
    cout << "     | || '_ \\ `--. \\/ _ \\/ _` | __|" << endl;
    cout << "    _| || | | /\\__/ /  __/ (_| | |_ " << endl;
    cout << "    \\___/_| |_\\____/ \\___|\\__, |\\__|" << endl;
    cout << "                             | |    " << endl;
    cout << "                             |_|    " << endl;
  }
  
  void ParserCL(int argc, char** argv) {
    int i;
    m_exe = argv[0];
    //cout << "Running " << m_exe << endl;
    for (i=strlen(argv[0])-1; i>=0; i--) {
      if (argv[0][i] == '/') {
	char tmp[1024];
	strcpy(tmp, argv[0]);
	tmp[i+1] = 0;
	m_exe = tmp;
      }
    }
    //cout << m_exe << endl;
    for (i=0; i<m_program.isize(); i++)
      m_program[i] = m_exe + m_program[i];
  
    bool bFull = false;
    if (argc == 2 && strcmp(argv[1], "-full") == 0)
      bFull =  true;
    
    bool bHelp = false;
    if (argc == 2 && strcmp(argv[1], "-h") == 0)
      bHelp =  true;
    if (argc == 2 && strcmp(argv[1], "-help") == 0)
      bHelp =  true;
    
    if (bFull || argc == 1 || bHelp) {
      if (!bHelp)
	cout << "Use the -h or -help option to print this screen." << endl << endl;
      ShowHelp(bFull);
      return;
    }
    int index = Find(argv[1]);
    if (index < 0) {
      cout << "Do not understand command '" << argv[1] << "'" << endl;
      cout << "Use the -h or -help option to show all available commands." << endl;
      return;
    }

    string doIt = m_program[index];
    for (i=2; i<argc; i++) {
      doIt += " ";
      doIt += argv[i];
    }
    //cout << "Run " << doIt << endl;
    int r = system(doIt.c_str());
  }
    

  void ShowHelp(bool bFull) {
    PrintLogo();
    cout << endl;
    int i;
    if (!bFull)
      cout << "For detailed options, type InSeqt -full" << endl << endl;
    cout << "Printing a list of available commands." << endl << endl;
    for (i=0; i<m_cmd.isize(); i++) {
      if (bFull)
	cout << "--------------------------------------------------";
      cout << endl << "Usage: InSeqt " << m_cmd[i] << " <options>" << endl;
      string cc = m_program[i] + " | grep description";
      if (bFull)
	cc = m_program[i];
      //cout << "Run " << cc << endl;
      int r = system(cc.c_str());
      //cout << endl;
    }
    cout << endl << endl;
  }

private:
  string m_exe;
  string m_args;

  int Find(const string & cmd) {
    for (int i=0; i<m_cmd.isize(); i++) {
      if (m_cmd[i] == cmd)
	return i;
    }
    return -1;
  }

  svec<string> m_cmd;
  svec<string> m_program;
  svec<string> m_depend;
};



int main( int argc, char** argv )
{
  Command cmd;
  cmd.ParserCL(argc, argv);

  return 0;
}
