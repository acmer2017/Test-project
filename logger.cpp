// Copyright 2017 Alan Kuhnle.

// This file is part of mim.

// mim is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// mim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with mim.  If not, see <http://www.gnu.org/licenses/>.
#ifndef logger_cpp
#define logger_cpp
#include <fstream>
#include <string>
#include <ctime>

using namespace std;

enum LogType { ERROR = -3, WARN = -2, OUTPUT = -1, INFO = 0, DEBUG = 1, TRACE = 2 };

class Logger {
public:
  LogType loglevel;
  ostream& of;
  bool echo;
  string msg;
  LogType msglevel;
  
  Logger() : of (cout) {
    //    of.open("log.txt");
    //    of = cout;
    loglevel = INFO;
    msglevel = INFO;
    echo = false;
  }
  
   Logger( LogType inlevel, ostream& os, bool echo_in = false ): of( os ) {
    //of.open("log.txt");
    //    of = cout;
    loglevel = inlevel;
    msglevel = INFO;
    echo = echo_in;
  }

  void set_level( LogType Lin ) {
    loglevel = Lin;
  }

  void operator()( LogType level, string msg ) {
    if (level <= loglevel) {
      switch (level) {
      case ERROR:
	of << time(NULL) << "\033[31m [ERROR] \033[0m" << msg << endl;
	break;
      case WARN:
	of << time(NULL) << "\033[33m [WARN] \033[0m" << msg << endl;
	break;
      case OUTPUT:
	of << time(NULL) << " [OUTPUT] " << msg << endl;
	break;
      case INFO:
	of << time(NULL) << "\033[32m [INFO] \033[0m" << msg << endl;
	break;
      case DEBUG:
	of << time(NULL) << "\033[31m [DEBUG] \033[0m" << msg << endl;
	break;
      case TRACE:
	of << time(NULL) << "\033[31m [TRACE] \033[0m" << msg << endl;
	break;
      }
      if (echo) {
	 if (level <= loglevel) {
	    switch (level) {
	    case ERROR:
	       cout << time(NULL) << "\033[31m [ERROR] \033[0m" << msg << endl;
	       break;
	    case WARN:
	       cout << time(NULL) << "\033[33m [WARN] \033[0m" << msg << endl;
	       break;
	    case OUTPUT:
	       cout << time(NULL) << " [OUTPUT] " << msg << endl;
	       break;
	    case INFO:
	       cout << time(NULL) << "\033[32m [INFO] \033[0m" << msg << endl;
	       break;
	    case DEBUG:
	       cout << time(NULL) << "\033[31m [DEBUG] \033[0m" << msg << endl;
	       break;
	    case TRACE:
	       cout << time(NULL) << "\033[31m [TRACE] \033[0m" << msg << endl;
	       break;
	    }
	 }
      }
    }	     
  }

  ~Logger() {
    //of.close();
  }
};

class endlclass {
  
} endL;

template <typename T>
Logger& operator<<( Logger& lhs, T word ) {
  lhs.msg += to_string(word);
  return lhs;
}

Logger& operator<<( Logger& lhs, const char* word ) {
  string tmp( word );
  lhs.msg += tmp;
  return lhs;
}

Logger& operator<<( Logger& lhs, string word ) {
  lhs.msg += word;
  return lhs;
}

Logger& operator<<( Logger& lhs, endlclass word ) {
  lhs( lhs.msglevel, lhs.msg );
  lhs.msg.clear();

  return lhs;
}

Logger& operator<<( Logger& lhs, LogType inmsgLevel ) {
  lhs.msglevel = inmsgLevel;

  return lhs;
}



#endif
