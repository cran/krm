#ifndef U_H
#define U_H

// added at the suggestion of Prof. Ripley
#include <stdio.h>

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>

#include <R.h>

using namespace std;

#define _PI 3.1415926535897932384626433832795
extern "C" { double myunif_rand(); }
//double myunif_rand(); 

void Tokenize(const string& str, vector<string>& tokens, const string& delimiters = " ");

class U
{

public:
	
	template <class T>
	static void copyArray(const T *from, T *to, int size) {
		for (int i=0; i<size; i++) to[i]=from[i];
	}

	template <class T>
	static bool compareArray(const vector<T>& a1, const vector<T>& a2, int size) {
		bool sameness = a1[0] == a2[0];
		bool matching=true;
		for (int i=1; i<size; i++) 
			if ( (a1[i]==a2[i]) != sameness) {
				matching=false;
				break;
			}
			
		return matching;
	}

	static void minmax(double a[], int size, double* minVal, double* maxVal) {
	    double _maxVal = a[0]; 
	    double _minVal = a[0]; 
	    for (int i=1; i<size; i++) {
	        if (a[i] > _maxVal) _maxVal = a[i];
	        if (a[i] < _minVal) _minVal = a[i];
		}
		*minVal=_minVal;
		*maxVal=_maxVal;
	}	

	template <class T>
	static T arraymax(T a[], int size) { // if named max, it is confused with system function max for some reason
	    T maxVal = a[0];     
	    for (int i=1; i<size; i++) 
	        if (a[i] > maxVal) maxVal = a[i];
	    return maxVal;
	}	

	template <class T>
	static T arraymax(std::vector<T> & a) { // if named max, it is confused with system function max for some reason
	    T maxVal = a[0];     
	    int size = a.size();
	    for (int i=1; i<size; i++) 
	        if (a[i] > maxVal) maxVal = a[i];
	    return maxVal;
	}	

	template <class T>
	static T arraymin(T a[], int size) {
	    T minVal = a[0];     
	    for (int i=1; i<size; i++) 
	        if (a[i] < minVal) minVal = a[i];
	    return minVal;
	}	

	template <class T>
	static T sum(T a[], int size){
	    T out = a[0];     
	    for (int i=1; i<size; i++) out += a[i];
	    return out;
	}

	template <class T>
	static T prod(T a[], int size){
	    T out = a[0];     
	    for (int i=1; i<size; i++) out *= a[i];
	    return out;
	}

	static double factorial(double a){		
	    return a>1?a*factorial(a-1):1;
	}

	template <class T>
    static string n2s ( T number ){
	  ostringstream oss;
	  oss<< number;
	  return oss.str();
	}

	static bool endsWith(string longstr, string str) {
		if (longstr.substr(longstr.length()-str.length(), str.length())==str) return true;
		else return false;
	}

	static void openWrite (ofstream& file, string fileName, ostream& logFile, bool overwrite=false) {
		if (!overwrite) {
			// should heck if file exists or not
			ifstream ifile(fileName.c_str());
			if (!ifile.fail()) {
				//The file exists, return error
				//if (logFile!=cout) cout << "File already exists: " << fileName << endl; // comment out b/c R CMD check complains
				logFile << "File already exists: " << fileName << endl;
				//exit (-1);			  // comment out b/c R CMD check complains
			} else {
				ifile.close();
			}
		}
		file.open(fileName.c_str(), ios::out);
		if(file.fail()) {
			//if (logFile!=cout) cout << "Fail to write " << fileName << endl; // comment out b/c R CMD check complains
			logFile << "Fail to write " << fileName << endl;
			//exit (-1);// comment out b/c R CMD check complains
		}
		else
			logFile << "  write " << fileName << endl;
	}

	static void openRead (ifstream& file, string fileName, ostream& logFile) {
		file.open(fileName.c_str(), ios::in);
		if(file.fail())
		{
			logFile << "Fail to read " << fileName << endl;
			//if (logFile!=cout) cout << "Fail to read " << fileName << endl;// comment out b/c R CMD check complains
			//exit (-1); // comment out b/c R CMD check complains
		}
		else
		{
			logFile << "  read " << fileName << endl;
		}
	}
	
	// read data from a matrix file, allocate memory for it on the heap
	template <class T>
	static void readData(string fileName, int& rowCount, int& columnCount, T**& data, ostream& logFile){		
		ifstream file1;
		U::openRead(file1, fileName, logFile);
		
		// find out how many columns there are
		char line_cstr[5000];
		file1.getline(line_cstr, 5000);
		string line(line_cstr);		
		vector<string> tokens;
    	Tokenize(line, tokens);
    	columnCount=tokens.size();
    	
		// find out how many rows there are
    	rowCount=0; 
    	while (!file1.eof()) {
			// if not end of file, increase by 1, this line should not come after file1.getline call
			rowCount++;
			file1.getline(line_cstr, 5000);
		}
		
		// close file, open again, and read data
		file1.close();
		U::openRead(file1, fileName, logFile);
		    	
	    data = new T*[rowCount];
	    for (int i=0; i<rowCount; i++) 
	        data[i] = new T[columnCount];
	
	    for (int i=0; i<rowCount; i++) {
	        for (int t=0; t<columnCount; t++)
	            file1 >> data[i][t];
	    }
	    file1.close();
		
		//U::print2Dmatrix(data, rowCount, columnCount, logFile);
	}

	// in the file, first line is header, the rest is a data matrix
	static void readData(string fileName, string* columnNames, int columnCount, long rowCount, double**out){
		ifstream file1;
        std::ostringstream oss;
		U::openRead(file1, fileName, oss); // change from cout to oss b/c R CMD check complains
		char line_cstr[500];
		file1.getline(line_cstr, 500);
		string line(line_cstr);
		
		vector<string> tokens;
    	Tokenize(line, tokens);
    	int totalColumnCount=tokens.size();
		
		std::vector<int> columnIdx(columnCount);
		for (int i=0; i<columnCount; i++) {
			string columnName=columnNames[i];
			for (int j=0; j<(int)(tokens.size()); j++) {
				if (columnName.compare(tokens[j])==0) {
					columnIdx[i]=j;
					break;
				}
				if (j==(int)(tokens.size()-1)) {
					//cout << columnName << " not found in header\n"; // comment out b/c R CMD check complains
					//exit(-1);// comment out b/c R CMD check complains
				}
			}
		}
		//U::printArray(columnIdx, columnCount);
		
		std::vector<double> numbers(totalColumnCount);
		for (long i=0; i<rowCount; i++) {
			for (int j=0; j<totalColumnCount; j++) file1 >> numbers[j];
			for (int j=0; j<columnCount; j++) out[i][j] = numbers[columnIdx[j]];
		}
		//U::print2Dmatrix(out, rowCount, columnCount);
	}

	static int getdir (string dir, vector<string> &files, string startsWith="")
	{
	    DIR *dp;
	    struct dirent *dirp;
	    if((dp  = opendir(dir.c_str())) == NULL) {
	        //cout << "Error opening " << dir << endl;// comment out b/c R CMD check complains
	        //exit(-1);// comment out b/c R CMD check complains
	    }
	
		string fileName;
	    while ((dirp = readdir(dp)) != NULL) {
			fileName=string(dirp->d_name);
			if (fileName.substr(0, startsWith.size())==startsWith) files.push_back(fileName);
	    }
	    closedir(dp);
	    return 0;
	}

	static void sendMailOnUnix(string to, string body) {
		char cwd[100];
		getcwd(cwd, 100);
		string scwd(cwd);
		
		char hostname[100];
		FILE *fhostname = popen("hostname", "r");
		fgets(hostname, 100, fhostname);
		string shostname(hostname);
		
		string subj=shostname.substr(0, shostname.length()-1)+scwd;
		
		string cmd="/usr/bin/mail -s '"+subj+"' "+to;
		FILE *mailer = popen(cmd.c_str(), "w");
		fprintf(mailer, "%s", body.c_str());
		pclose(mailer); 
	}

	static int finishUp (int argc, char *argv[], string _msg, time_t& t0, ostream& logFile, string email="") {
		stringstream msg;
		msg << endl;
		
		int minPassed = ((long)time (NULL)-t0)/60;
		int dayPassed = minPassed/60/24;
		int hourPassed = minPassed/60-dayPassed*24;
		for (int i=0; i<argc; i++) msg<< argv[i] << " ";
		msg << endl;
		
		struct tm *current = localtime(&t0);
		msg << "Started at 200" << current->tm_year-100 <<"/"<< current->tm_mon+1 <<"/"<< current->tm_mday <<" "
		    << current->tm_hour <<":"<< current->tm_min <<":"<< current->tm_sec << endl;

		msg << "Finished in ";
		if (dayPassed>0) msg << U::n2s(dayPassed) << " days, ";
		if (dayPassed>0 ||hourPassed>0) msg << U::n2s(hourPassed) << " hours, ";
		msg << U::n2s(minPassed%60)+" min and " << ((long)time (NULL)-t0)%60 << " sec, or "
			<< U::n2s((long)time (NULL)-t0) << " seconds.\n";
		
		msg << _msg << endl;

		logFile << msg.str(); 
		//if (logFile!=cout) cout << msg.str();// comment out b/c R CMD check complains

		if (email!="" and minPassed > 5) 
			U::sendMailOnUnix(email, msg.str());
		
		#ifdef __DEBUG__		
		//cout << "\nDone. Please press ENTER to continue..." << endl;
		//cin.get();
		#endif		
		return 1;
	}
	
	template <class T>
	static void printArray(const T* array, int dim, ostream& save, bool println=true) {
		for (int i=0; i<dim; i++) save << array[i] << " ";
		if (println) save << endl;
	}

	template <class T>
	static void print2Dmatrix(T** matrix, int dim1, int dim2, ostream& save) {
		for (int i=0; i<dim1; i++) {
			for (int j=0; j<dim2; j++) save << matrix[i][j] << " ";
			save << endl;
		}
	}

	template <class T>
	static void print2Dmatrix(T* matrix, int dim1, int dim2, ostream& save) {
		for (int i=0; i<dim1; i++) {
			for (int j=0; j<dim2; j++) save << matrix[i*dim2+j] << " ";
			save << endl;
		}
	}
	
	template <class T>
	static void printArray (vector<T> v, ostream& save) {
		for (int i=0; i<v.size(); i++) save << v[i] << " ";
		save << endl;
	}
	
	// return the log of the sum of weight*exp(logValues)
	// todo: need to deal with max and min in a weight-sensitive way
	static double logSumExp (double * logValues, int howMany, double* weights){
		if (howMany==1) return logValues[0]+log(weights[0]);

		double maximum=arraymax(logValues, howMany);
		if (maximum==R_NegInf) return R_NegInf;
		if (maximum==R_PosInf) return R_PosInf;

		double *shiftedLogValues = new double [howMany];
		for (int i=0; i<howMany; i++) shiftedLogValues[i]=logValues[i]-maximum;
		double result=0;		
		for (int i=0; i<howMany; i++) result += exp(shiftedLogValues[i]) * weights[i]; 
		delete[] shiftedLogValues;
		return log(result)+maximum;
	}
	static double logSumExp (std::vector<double> & logValues, double* weights){
           int howMany = logValues.size();
		if (howMany==1) return logValues[0]+log(weights[0]);

		double maximum=arraymax(logValues);
		if (maximum==R_NegInf) return R_NegInf;
		if (maximum==R_PosInf) return R_PosInf;

		double *shiftedLogValues = new double [howMany];
		for (int i=0; i<howMany; i++) shiftedLogValues[i]=logValues[i]-maximum;
		double result=0;		
		for (int i=0; i<howMany; i++) result += exp(shiftedLogValues[i]) * weights[i]; 
		delete[] shiftedLogValues;
		return log(result)+maximum;
	}
	// todo: need to deal with max and min in a weight-sensitive way
	static double logSumExpLogW (double* logValues, int howMany, double* logWeights){
		if (howMany==1) return logValues[0]+logWeights[0];

		double maximum=arraymax(logValues, howMany);
		if (maximum==R_NegInf) return R_NegInf;
		if (maximum==R_PosInf) return R_PosInf;

		double *shiftedLogValues = new double [howMany];
		for (int i=0; i<howMany; i++) shiftedLogValues[i]=logValues[i]-maximum;
		double result=0;		
		for (int i=0; i<howMany; i++) result += exp(shiftedLogValues[i] + logWeights[i]); 
		delete[] shiftedLogValues;
		return log(result)+maximum;
	}
	static double logSumExp (double* logValues, int howMany){
		if (howMany==1) return logValues[0];
		double *weights = new double [howMany];
		for (int i=0; i<howMany; i++) weights[i]=1;
		double result = logSumExp(logValues, howMany, weights);
		delete[] weights;
		return result;
	}
	static double logAddExp (double logV1, double logV2){
		double logValues[2];
		logValues[0]=logV1;
		logValues[1]=logV2;
		return logSumExp(logValues, 2);
	}
	// return log ( exp(logV1) - exp(logV2) ), but make sure logV1 is larger than logV2
	static double logSubstractExp (double logV1, double logV2){
		return -1;
		// to be done
	}
	
	// z is a vector of group indicators from 0 to k-1 of length n
	// count is a vector of counts of length k
	// return number of non-empty components
	template <class T>
	static int table(T* z, int n, int * count, int k) {
		for (int j=0; j<k; j++) count[j]=0;
		for (int i=0; i<n; i++) count[ z[i] ] ++;
		int truek=0;
		for (int j=0; j<k; j++) 
			if (count[j]!=0) truek++;
		return truek;
	}
	template <class T>
	static int table(vector<T>& z, int * count, int k) {
		int n=z.size();
		for (int j=0; j<k; j++) count[j]=0;
		for (int i=0; i<n; i++) count[ z[i] ] ++;
		int truek=0;
		for (int j=0; j<k; j++) 
			if (count[j]!=0) truek++;
		return truek;
	}
	
	// map a string to an integer, could be useful when seeding a rng
	static unsigned int stringHashCode (string s) {
		unsigned int out=0;
		for (size_t i=0; i<s.size(); i++)
			out+=(unsigned int)s[i];
		return out;
	} 

};

# endif
