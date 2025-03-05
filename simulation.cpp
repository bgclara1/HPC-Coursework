#include <iostream>
#include <map>
#include <unordered_map>

using namespace std;

int main(int argc, char *argv[]) {              //read cmd args w main params.
    int i = 0;

    // test reading in args properly
    while (i < argc) {  //argc num args provided
        cout << "Argument " << i + 1 << ": " << argv[i]
             << endl;
        i++;
    }
    
    map<char, double> testCaseDict;
   // testCaseDict['-ic-one'] = [1.0,10,10,10,0.0,0.0,0.0,0];



    return 0;
}