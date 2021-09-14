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
#ifndef MCKP_CPP
#define MCKP_CPP

//input is system of profits, each
//row must have the same length
// 
//the costs are assumed to be 
//a( x_{ij} ) = j, since for influence
//max, the cost is number of seed nodes,
//for our application

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

namespace mckp {

   class mypair {
   public: 
      vector< unsigned > x;
      double c_x;

      mypair& operator= ( mypair const &other ) {
	 c_x = other.c_x;
	 x.clear();
	 x.assign( other.x.begin(), other.x.end() );

	 return *this;
      }

      mypair() {
	 c_x = 0;
      }

      mypair( const mypair& other ) {
	 c_x = other.c_x;
	 x.clear();
	 x.assign( other.x.begin(), other.x.end() );
      }
   };

   bool operator== (const mypair& lhs, const mypair& rhs ) {

      return (lhs.c_x == rhs.c_x);
      // for (unsigned i = 0; i < rhs.x.size(); ++i) {
      //   if (i >= lhs.x.size()) {
      //     return false;
      //   }
      //   if ((lhs.x)[i] != (rhs.x)[i]) {
      //     return false;
      //   }
      // }

   }

   bool operator< (const mypair& lhs, const mypair& rhs) {
      return (lhs.c_x < rhs.c_x);
   }

   bool operator<= (const mypair& lhs, const mypair& rhs) {
      return (operator< (lhs, rhs) || operator==(lhs, rhs));
   }

   class item {
   public:
      double c;
      double a;
      unsigned i;
      unsigned j;
      bool x; //included in system

      item& operator=( item const &other ) {
	 this -> c = other.c;
	 this -> a = other.a;
	 this -> i = other.i;
	 this -> j = other.j;
	 this -> x = other.x;

	 return *this;
      }
  
      item() {
	 this -> c = 0;
	 this -> a = 0;
	 this -> i = 0;
	 this -> j = 0;
	 this -> x = 0;
      }

      item( const item& other ) {
	 this -> c = other.c;
	 this -> a = other.a;
	 this -> i = other.i;
	 this -> j = other.j;
	 this -> x = other.x;
      }
   };

   bool operator< ( const item& lhs, const item& rhs) {
      return (lhs.c / lhs.a) < (rhs.c / rhs.a);
   }
   bool operator<= (const item& lhs, const item& rhs) {
      return (lhs.c / lhs.a) <= (rhs.c / rhs.a);
   }
   bool operator== (const item& lhs, const item& rhs) {
      bool val = (lhs.a == rhs.a) && (lhs.c == rhs.c) && (lhs.i == rhs.i) && (lhs.j == rhs.j);
      return val;
   }

   void transform( vector < vector < item > >& sys );
   void trim_sys( vector< vector< item > >& sys );
   void enforce_31( vector < vector < item > >& sys );
   void enforce_32( vector < vector < item > >& sys );

   void print_sys( vector < vector < item > >& sys, ostream& of ) {
      of << "System:\n";
      for (unsigned i = 0; i < sys.size(); ++i) {
	 for (unsigned j = 0; j < sys[i].size(); ++j) {
	    of << "(" << sys[i][j].c << ',' << sys[i][j].a << ')' << "  ";
	 }
	 of << endl;
      }
   }

   void print_choice( vector< unsigned >& choice, ostream& os ) {
      os << "Choice:\n";
      for (unsigned i = 0; i < choice.size(); ++i) {
	 os << i + 1 << ' ' << choice[i] << endl;
      }
   }

   void BM( vector<unsigned>& jprime, vector< vector< item > >& in_sys, double b ) {
      vector< vector< item > > sys( in_sys.begin(), in_sys.end() );
      //enforce Lemmas 3.1, 3.2, to get R'
      enforce_31( sys );
      trim_sys( sys );
      enforce_32( sys );
      trim_sys( sys );
      cout << "After lemmas:\n";
      print_sys( sys, cout );
      //transform R' to R''
      vector< vector< item > > osys( sys.begin(), sys.end() ); //save copy before transformation to R''
      transform( sys );
      cout << "After transform:\n";
      print_sys( sys, cout );

      // //map entries in sys to those in in_sys, so we can return intelligible result
      // vector< vector < unsigned > > mapping;
      // for (unsigned i = 0; i < sys.size(); ++i) {
      //   vector <unsigned> match_row;
      //   for (unsigned k = 0; k < sys[i].size(); ++k) {
      //     //find match for sys[i][k]
      //     for (unsigned j = 0; j < in_sys[i].size(); ++j) {
      //       if (sys[i][k].a == in_sys[i][j].a) {
      //         if (sys[i][k].c == in_sys[i][j].c) {
      //           //found match
      //           match_row.push_back( j );
      //           break;
      //         }
      //       }
      //     }
      //   }
      //   mapping.push_back( match_row );
      // }
      vector< item > copy_sys;
      unsigned k = sys.size(); //number of rows
      for (unsigned i = 0; i < k; ++i) {
	 for (unsigned j = 0; j < sys[i].size(); ++j) {
	    copy_sys.push_back( sys[i][j] );
	 }
      }

      sort( copy_sys.begin(), copy_sys.end() ); //sorts ascending order
      reverse( copy_sys.begin(), copy_sys.end() ); //we need descending

      unsigned r = 0;
      double sum = 0; //copy_sys[r].a;
      bool debug = true;
      while (sum <= b) {

	 ++r;    
	 sum += copy_sys[r - 1].a;

	 if (debug) 
	    cerr << sum << ' ' << r << ' ' << '(' << copy_sys[r - 1].c << ',' 
		 << copy_sys[r - 1].a << ')' << endl;

      }
      r -= 1;

      if (debug) {
	 sum = 0;
	 for (unsigned kk = 0; kk < r; ++kk) {
	    sum += copy_sys[kk].a;
	 }
	 cout << sum << endl;
	 cout << "r: " << r << endl;
	 cout << '(' << copy_sys[r].i << ',' << copy_sys[r].j << ')' << endl;
	 for (unsigned i = 0; i < copy_sys.size(); ++i) {
	    cerr << '(' << copy_sys[i].i << ',' 
		 << copy_sys[i].j << ')' << ' ';
	 }
	 cerr << endl;

	 for (unsigned i = 0; i < copy_sys.size(); ++i) {
	    cerr << '(' << copy_sys[i].c << ',' 
		 << copy_sys[i].a << ')' << ' ';
	 }
	 cerr << endl;
      }
  
      jprime.clear();
      jprime.assign( k , 0 ); 
      //each row gets a jprime
  
      for (unsigned i = 0; i < k; ++i) { 
	 //iterate through each row
	 for (unsigned j = 0; j < r; ++j) {
	    //iterate through first r pairs, starting at 0
	    if (copy_sys[j].i == (i + 1)) { //row label is +1
	       if (copy_sys[j].j > jprime[i]) {
		  jprime[i] = copy_sys[j].j;
	       }
	    }
	 }
      }

      print_choice( jprime, cout );
  
      double pval = 0;
      for (unsigned i = 0; i < k; ++i) {
	 if (jprime[i] != 0) {
	    pval += osys[i][ jprime[i] - 1 ].c;
	 }
      }
      cout << "Pval\n";
      cout << pval << endl;
      double alt_pval = in_sys[copy_sys[r].i - 1][copy_sys[r].j -1 ].c;
      cout << alt_pval << endl;
  
      if (pval < alt_pval ) {
	 jprime.assign(k, 0);
	 jprime[ copy_sys[r].i - 1 ] = copy_sys[r].j;
      }

      //map choice back to original system
  
      for (unsigned i = 0; i < jprime.size(); ++i) {
	 if (jprime[i] != 0) 
	    jprime[i] = osys[i][jprime[i] - 1].a;
      }

      if (debug) {
	 // print_sys( in_sys, cout );
	 print_choice( jprime, cout );
      }
   }

   bool next_i( vector< vector< unsigned > >& result, 
		unsigned n, unsigned m ) {
      result.clear();
      if (n == 0) {
	 return false;
      }
      if (n == 1) {
	 //base case
	 vector< unsigned > z(1, 0);
	 result.push_back( z );
	 if (m == 0) {
      
	 } else {
	    vector< unsigned > z(1, 1);
	    result.push_back( z );
	 }

	 return true;
      }

      if (m == 0) {
	 vector< unsigned > z( n, 0 );
	 result.push_back( z );
	 return true;
      }
      if (m >= n)
	 m = n;
    
      //recursive case
      vector< vector< unsigned > > smaller_res;
      vector< unsigned > yright; 
      vector< unsigned > yleft;
      next_i( smaller_res, n - 1, m - 1 );

      unsigned i = 0;
      //  for (unsigned i = 0; i < n; ++i ) {
      for (unsigned j = 0; j < smaller_res.size(); ++j) {
	 vector< unsigned > x = smaller_res[j];
	 yleft.clear();
	 yleft.push_back( 1 );
      
	 vector< unsigned > z;
	 z.insert( z.end(), yleft.begin(), yleft.end() );
	 z.insert( z.end(), x.begin(),x.end() );
	 result.push_back( z );
	 if (z.size() > n) {
	    cerr << "Error: " <<  n << ' ' << m << ' ' << z.size() << endl;
	 }
	 if (z.size() < n) {
	    cerr << "Error: " <<  n << ' ' << m << ' ' << z.size() << endl;
	 }
      }
      //  }
  
      smaller_res.clear();
      next_i( smaller_res, n - 1, m  );

      for (unsigned j = 0; j < smaller_res.size(); ++j) {
	 vector< unsigned > x = smaller_res[j];
	 yleft.clear();
	 yleft.push_back( 0 );
      
	 vector< unsigned > z;
	 z.insert( z.end(), yleft.begin(), yleft.end() );
	 z.insert( z.end(), x.begin(), x.end() );
	 result.push_back( z );
      }
      //  }

      return true;

   }

   double cost( vector< unsigned>& x, 
		vector< vector< item > >& sys ) {
      unsigned k = sys.size();
      unsigned l = sys[0].size();

      double cost = 0.0;

      for (unsigned i = 0; i < k; ++i) {
	 for (unsigned j = 0; j < l; ++j) {
	    cost += sys[i][j].a * x[ i * l + j ];
	 }
      }

      return cost;
   }

   double profit( vector< unsigned>& x, 
		  vector< vector< item > >& sys ) {
      unsigned k = sys.size();
      unsigned l = sys[0].size();

      double profit = 0.0;

      for (unsigned i = 0; i < k; ++i) {
	 for (unsigned j = 0; j < l; ++j) {
	    profit += sys[i][j].c * x[ i * l + j ];
	 }
      }

      return profit;
   }

   void make_sysp(vector< vector < item > >& sys,
		  vector< unsigned >& x,
		  double b,
		  vector< vector < item > >& sysp,
		  double& bp ) {
      unsigned k = sys.size();
      unsigned l = sys[0].size();

      bp = b - cost(x, sys);
      double chat = profit(x, sys);
      vector < bool > flagrow (k, false);
      for (unsigned i = 0; i < k; ++i ) {
	 for (unsigned j = 0; j < l; ++j ) {
	    if (x[i*l + j] == 1) {
	       flagrow[i] = true;
	       if ( sys[i][j].c < chat ) {
		  chat = sys[i][j].c;
	       }
	    }
	 }
      }

      vector< item > new_row;
      sysp.clear();
      for (unsigned i = 0; i < sys.size(); ++i) {
	 new_row.clear();
	 for (unsigned j = 0; j < sys[i].size(); ++j ) {
	    item tmp = sys[i][j];
	    if (flagrow[i])
	       tmp.c = 0;
	    if (tmp.c > chat)
	       tmp.c = 0;
	    new_row.push_back( tmp );
	 }
	 sysp.push_back( new_row );
      }
  
   }
   bool feasible( vector< unsigned >& x, 
		  vector< vector< item > >& sys,
		  double b ) {
      if (cost( x, sys ) > b) {
	 return false;
      }
      unsigned k = sys.size();
      unsigned l = sys[0].size();
      double picked;
      for (unsigned i = 0; i < k; ++i) {
	 picked = 0.0;
	 for (unsigned j = 0; j < sys[i].size(); ++j ) {
	    picked += x[ l * i + j ];
	 }
	 if (picked > 1.0)
	    return false;
      }

      return true;
   }


   void BM_epsi( vector<unsigned>& choice, vector< vector< item > >& sys, double b, double epsilon ) {
      unsigned k = sys.size();
      unsigned l = sys[0].size();
      vector< vector < unsigned > > x;
      unsigned n = k * l;
      unsigned m = (unsigned) (1/epsilon - 1);
      cerr << "n, m: " << n << ',' <<  m << endl;
      next_i( x, n, m );
      bool debug = true;
      cout << "x\n";
      unsigned xcnt = 0;
      if (debug) {
	 for (unsigned i = 0; i < x.size(); ++i) {
	    for (unsigned j = 0; j < x[i].size(); ++j) {
	       cout << x[i][j];
	       ++xcnt;
	    }
	    cout << endl;
	 }
      }

      cout << xcnt << endl;
  
      vector< vector< item > > sysp;
      double bp;
      vector< mypair > psols;
      for (unsigned r = 0; r < x.size(); ++r) {
	 if ( feasible( x[r], sys, b ) ) {
	    make_sysp( sys, x[r], b, sysp, bp );
	    vector< unsigned > xp;
	    BM( xp, sysp, bp );
	    //need to convert xp to a vector of length x[r]
	    //xp is a vector of length k
	    //such that if no item is picked from row i
	    //xp[i] = 0 and xp[i] = j + 1 if item j picked
	    vector< unsigned > xpp( n, 0 );
	    cout << "mark1" << endl;
	    cout << n << ' ' << k << ' ' << xp.size() << endl;
	    cout << l << endl;
	    for (unsigned ii = 0; ii < k; ++ii) {
	       if (xp[ii] != 0) {
		  cout << xp[ii] << endl;
		  cout << ii*l + (xp[ii] - 1) << endl;
		  xpp[ ii * l + (xp[ii] - 1) ] = 1;
	       }
	    }
	    cout << "profit1" << endl;
	    double c_1 = profit(x[r], sys);
	    cout << "profit2" << endl;
	    double c_2 = profit(xpp, sysp);
	    mypair tmp;
	    tmp.x.assign( x[r].begin(), x[r].end() );
	    if (c_1 >= c_2) 
	       tmp.c_x = c_1;
	    else
	       tmp.c_x = c_2;
	    psols.push_back( tmp );

	    for (unsigned ii = 0; ii < x[r].size(); ++ii) {
	       cout << x[r][ii];
	    }
	    cout << endl;
	    cout << c_1 << endl;
	    for (unsigned ii = 0; ii < xp.size(); ++ii) {
	       cout << xp[ii];
	    }
	    cout << endl;
	    for (unsigned ii = 0; ii < xpp.size(); ++ii) {
	       cout << xpp[ii];
	    }
	    cout << endl << c_2 << endl << endl;
	 }
      }

      mypair sol = psols[0];
      for (unsigned i = 0; i < psols.size(); ++i) {
	 if (psols[i].c_x > sol.c_x) {
	    sol = psols[i];
	 }
      }

      choice.assign( sol.x.begin(), sol.x.end() );
   }

   void read_input( istream& ifs, vector < vector< item > >& sys, double& b ) {
      //  ifstream ifs( ifile.c_str() );
      string s;
      int i = -1;
      stringstream ss;
      while (getline( ifs, s ) ) {
	 ss.clear();
	 ss.str( s );
	 if (i == -1) {
	    //first line is b
	    ss >> b;
	 } else {

	    unsigned j = 0;    
	    double tmp;
	    vector < item> row;
	    while( ss >> tmp ) {

	       item i_tmp;
	       i_tmp.c = tmp;
	       i_tmp.a = j + 1;
	       i_tmp.i = i + 1;
	       i_tmp.j = j + 1;
	       i_tmp.x = true;

	       row.push_back( i_tmp );
	       ++j;
	    }

	    sys.push_back( row );
	 }
	 ++i;
      }

   }

   void enforce_31( vector < vector < item > >& sys ) {
      bool debug = true;
      if (debug) {
	 cout << "enforce_31\n";
	 print_sys( sys, cout );
      }
      //need to check each row to see
      //if we may remove an element
      //which is less profit effective
      //and less profit than another
      unsigned k = sys.size();
      for (unsigned s = 0; s < k; ++s) {
	 //check each pair on this row
	 for (unsigned p = 0; p < sys[s].size(); ++p) {
	    if (sys[s][p].x) {
	       for (unsigned q = p + 1; q < sys[s].size(); 
		    ++q) {
		  if (sys[s][q].x && sys[s][q].a > 0) {
		     if (sys[s][p].x && sys[s][p].a > 0 && sys[s][p].c > 0) {
			//determine if we need to discard one
			if (sys[s][p].c >= sys[s][q].c) {
			   double t1 = sys[s][p].c / sys[s][p].a;
			   double t2 = sys[s][q].c / sys[s][q].a;
			   if (t1 >= t2) {
			      //we can discard q
			      sys[s][q].x = false;
			   }
			} else {
			   double t1 = sys[s][p].c / sys[s][p].a;
			   double t2 = sys[s][q].c / sys[s][q].a;
			   if (t2 >= t1) {
			      //we can discard p
			      sys[s][p].x = false;
			   }
			}
		     }
		  }
	       }
	    }
	 }
    
      }

      if (debug) {
	 print_sys( sys, cout );
	 cout << "end enforce_31" << endl;
      }
  
   }

   void enforce_32( vector < vector < item > >& sys ) {
      bool debug = true;
      if (debug) {
	 cout << "enforce_32" << endl;
	 print_sys( sys, cout );

      }
      //need to check each row to see
      //if we may remove an element
      //which is less profit effective
      //and less profit than another
      unsigned k = sys.size();
      for (unsigned s = 0; s < k; ++s) {
	 //check each pair on this row
	 for (unsigned p = 0; p < sys[s].size(); ++p) {
	    if (sys[s][p].x) {
	       for (unsigned r = p + 2; r < sys[s].size(); ++r) {
		  for (unsigned q = p + 1; q < r; ++q ) {
		     if (sys[s][q].x) {
			if (sys[s][p].x) {
			   if (sys[s][r].x) {
			      //determine if we need to discard one
			      if (sys[s][q].c == sys[s][p].c) {
				 sys[s][q].x = false;
			      } else {
				 if (sys[s][r].a == sys[s][q].a) {
				    sys[s][q].x = false;
				 }
				 else {
				    double t1 = (sys[s][q].c - sys[s][p].c) / (sys[s][q].a - sys[s][p].a);
				    double t2 = (sys[s][r].c - sys[s][q].c) / (sys[s][r].a - sys[s][q].a);
				    if (t1 <= t2) {
				       sys[s][q].x = false;
				    }
				 }
			      }

                  
			   }

			}
		     }
		  }
	       }
	    }
	 }
      }

      if (debug) {
	 print_sys( sys, cout );
	 cout << "end enforce_32" << endl;
      }
  
   }



   void trim_sys( vector< vector< item > >& sys ) {
      vector< vector< item > > new_sys;
      unsigned new_j;
      for (unsigned i = 0; i < sys.size(); ++i) {
	 vector< item > new_row;
	 new_j = 1;
	 for (unsigned j = 0; j < sys[i].size(); ++j) {
	    if (sys[i][j].x) {
	       sys[i][j].j = new_j;
	       new_row.push_back( sys[i][j] );
	       ++new_j;
	    }
	 }
	 new_sys.push_back( new_row );
      }

      sys.clear();
      sys.assign( new_sys.begin(), new_sys.end() );
   }

   void transform( vector < vector < item > >& sys ) {
      vector< vector< item > > new_sys( sys.begin(), sys.end() );
      for (unsigned i = 0; i < sys.size(); ++i) {
	 for (unsigned j = 0; j < sys[i].size(); ++j) {
	    if (j != 0) {
	       new_sys[i][j].a = sys[i][j].a - sys[i][j - 1].a;
	       new_sys[i][j].c = sys[i][j].c - sys[i][j - 1].c;
	    }
	 }
      }

      sys.clear();
      sys.assign( new_sys.begin(), new_sys.end() );
   }


}

// int main(int argc, char** argv) {
//    bool debug = false;
//   //read input file
//   vector < vector< item > > sys;

//   double b;
//   string ifile( argv[1] );
//   read_input_file( ifile, sys, b );

//   //Run alg. BM
//   //iterate through costs b
//   cout << "b: " << b << endl;
//   vector< unsigned > choice;
//   //  BM(choice, sys, (double) b );
//   BM_epsi(choice, sys, (double) b, 0.1 );
//   //  cout << "Here?" << endl;
//   if(debug) {
//      for (unsigned i = 0; i < sys.size(); ++i) {
//         for (unsigned j = 0; j < sys[i].size(); ++j) {
// 	   cerr << sys[i][j].c << ' ';
//         }
//         cerr << endl;
//       }
//     }



//   //string ofname = "./out/out";
//   //    ofname = ofname + to_string(b) + ".txt";
//   //    ofstream ofile( ofname );

//   cout << "size of choice: " << choice.size() << endl;
//   for (unsigned i = 0; i < choice.size(); ++i) {
//        cout << i + 1 << ' ' << choice[i] << endl;
//        //     ofile << i + 1 << ' ' << choice[i] << endl;
//      }
//   //     ofile.close();
//   //     ofile.clear();

//   return 0;
   
// }


#endif
