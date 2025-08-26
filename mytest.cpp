#include "dnadb.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
using namespace std;
//Random class that was previously used                                                                                                                                           
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE};
class Random {
public:
  Random(){}
  Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
  {
    if (type == NORMAL){
      m_generator = std::mt19937(m_device());
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
    else if (type == UNIFORMINT) {
      m_generator = std::mt19937(10);// 10 is the fixed seed value                                                                                                                
      m_unidist = std::uniform_int_distribution<>(min,max);
    }
    else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers                                                                                            
      m_generator = std::mt19937(10);// 10 is the fixed seed value                                                                                                                
      m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
    }
    else {
      m_generator = std::mt19937(m_device());
    }
  }
  void setSeed(int seedNum){
    m_generator = std::mt19937(seedNum);
  }
  void init(int min, int max){
    m_min = min;
    m_max = max;
    m_type = UNIFORMINT;
    m_generator = std::mt19937(10);// 10 is the fixed seed value                                                                                                                  
    m_unidist = std::uniform_int_distribution<>(min,max);
  }
  void getShuffle(vector<int> & array){
    for (int i = m_min; i<=m_max; i++){
      array.push_back(i);
    }
    shuffle(array.begin(),array.end(),m_generator);
  }
  void getShuffle(int array[]){
    vector<int> temp;
    for (int i = m_min; i<=m_max; i++){
      temp.push_back(i);
    }
    std::shuffle(temp.begin(), temp.end(), m_generator);
    vector<int>::iterator it;
    int i = 0;
    for (it=temp.begin(); it != temp.end(); it++){
      array[i] = *it;
      i++;
    }
  }

  int getRandNum(){
    int result = 0;
    if(m_type == NORMAL){
      result = m_min - 1;
      while(result < m_min || result > m_max)
        result = m_normdist(m_generator);
    }
    else if (m_type == UNIFORMINT){
      result = m_unidist(m_generator);
    }
    return result;
  }
  double getRealRandNum(){
    double result = m_uniReal(m_generator);
    result = std::floor(result*100.0)/100.0;
    return result;
  }
  string getRandString(int size){
    string output = "";
    for (int i=0;i<size;i++){
      output = output + (char)getRandNum();
    }
    return output;
  }

  int getMin(){return m_min;}
  int getMax(){return m_max;}
private:
  int m_min;
  int m_max;
  RANDOM m_type;
  std::random_device m_device;
  std::mt19937 m_generator;
  std::normal_distribution<> m_normdist;
  std::uniform_int_distribution<> m_unidist;
  std::uniform_real_distribution<double> m_uniReal;
};
unsigned int simpleHash(string key) {
    unsigned int hash = 0;
    for (char ch : key) {
        hash = 29 * hash + ch;
    }
    return hash;
}

class Tester {
public:
    bool testInsertNoCollision();
    bool testInsertWithCollision();
    bool testRejectDuplicate();
    bool testInvalidLocation();
    bool testRemoveSuccess();
    bool testRemoveNonExistent();
    bool testRehashByLoadFactor();
    bool testRehashByDeleteRatio();
    bool testUpdateLocation();
    bool testPolicySwitch();
    bool testFindFromOldTable();
};
bool Tester::testInsertNoCollision() {//inserting two non-colliding DNA objects and getting them                                                                                  
    DnaDb db(101, simpleHash, LINEAR);
    DNA d1("ACGTA", 100000);
    DNA d2("TGCAA", 100001);

    if (!db.insert(d1)) {//if the first DNA object is not inserted, it shows lets you know it failed                                                                              
      cout << "Insert failed for d1" << endl;
        return false;
    }

    if (!db.insert(d2)) {//will do same here                                                                                                                                      
      cout << "Insert failed for d2" << endl;
        return false;
    }

    if (!(db.getDNA("ACGTA", 100000) == d1)) {//Tries not getting DNA for the first, will let you know it failed                                                                  
      cout << "GetDNA failed for d1" << endl;
        return false;
    }

    if (!(db.getDNA("TGCAA", 100001) == d2)) {//same here for the second DNA                                                                                                      
      cout << "GetDNA failed for d2" << endl;
        return false;
    }

    return true;
}

bool Tester::testInsertWithCollision() {
  //Makes a hash table of size 101 using simpleHash that is Linear probing                                                                                                        
  DnaDb db(101, simpleHash, LINEAR);
  //makes two dnas that are most likely going collide                                                                                                                             
  DNA d1("AAAAA", 100000);
  DNA d2("AAAAB", 100001);

  //insert both dnas even if they collide                                                                                                                                         
  if (!db.insert(d1)) return false;
  if (!db.insert(d2)) return false;
  //check that you can get both dnas                                                                                                                                              
  if (!(db.getDNA("AAAAA", 100000) == d1)) return false;
  if (!(db.getDNA("AAAAB", 100001) == d2)) return false;

  return true;
}

bool Tester::testRejectDuplicate() {
  DnaDb db(101, simpleHash, LINEAR);
  DNA d("ACGTA", 100000);
  if (!db.insert(d)) return false;//inserts dna                                                                                                                                   
  if (db.insert(d)) return false;//makes sure you cant enter duplicate                                                                                                            
  return true;
}

bool Tester::testInvalidLocation() {
  DnaDb db(101, simpleHash, LINEAR);
  DNA d("ACGTA", 99999);  // Invalid location ID                                                                                                                                  
  if (db.insert(d)) return false;
  return true;
}

bool Tester::testRemoveSuccess() {
  DnaDb db(101, simpleHash, LINEAR);
  DNA d("ACGTA", 100000);
  db.insert(d);
  if (!db.remove(d)) return false;//remove dna                                                                                                                                    
  if (!(db.getDNA("ACGTA", 100000) == d)) return false; // still returned by getDNA                                                                                               
  return true;
}

bool Tester::testRemoveNonExistent() {
    DnaDb db(101, simpleHash, LINEAR);
    DNA d("ACGTA", 100000);//create dna                                                                                                                                           
    if (db.remove(d)) return false;//try to remove it without inserting                                                                                                           
    return true;
}

bool Tester::testRehashByLoadFactor() {
  DnaDb db(101, simpleHash, LINEAR);

  //generate random for location and DNA chars                                                                                                                                    
  Random randId(100000, 999999, UNIFORMINT);
  Random randBase(0, 3, UNIFORMINT); //for DNA characters A/C/G/T                                                                                                                 
  string bases = "ACGT";

  vector<DNA> inserted;//storing DNA chars                                                                                                                                        

  for (int i = 0; i < 60; ++i) {//put more than allowed to force a rehash                                                                                                         
    string sequence = "";
    for (int j = 0; j < 5; ++j) sequence += bases[randBase.getRandNum()];
    int locId = randId.getRandNum();
    DNA d(sequence, locId);
    if (!db.insert(d)) {
      cout << "Insert failed at i = " << i << " with sequence = " << sequence << ", locId = " << locId << endl;
      return false;
    }
    inserted.push_back(d);
  }

  for (int i = 0; i < 30; ++i) {
    string sequence = "";
    for (int j = 0; j < 5; ++j) sequence += bases[randBase.getRandNum()];
    db.insert(DNA(sequence, randId.getRandNum()));
  }
  //checks to make sure all that were originally inserted are still accessible                                                                                                    
  for (int i = 0; i < (int)inserted.size(); ++i) {
    const DNA& d = inserted[i];
    DNA result = db.getDNA(d.getSequence(), d.getLocId());
    if (result.getSequence().empty()) {
      cout << "GetDNA failed for " << d.getSequence() << ", " << d.getLocId() << endl;
      return false;
    }
  }
  return true;
}

bool Tester::testRehashByDeleteRatio() {
  DnaDb db(101, simpleHash, LINEAR);
  vector<DNA> data;
  for (int i = 0; i < 50; ++i) {//put 50 dna at different locations                                                                                                               
    data.push_back(DNA("AACC" + to_string(i), 100000 + i));
    db.insert(data.back());
  }
  for (int i = 0; i < 45; ++i) {//remove 45 dna to force rehashing with >0.8                                                                                                      
    db.remove(data[i]);
  }
  for (int i = 45; i < 50; ++i) {//Check to see if you can still access the last 5                                                                                                
    if (db.getDNA(data[i].getSequence(), data[i].getLocId()).getSequence().empty()) {
      cout << "GetDNA failed at index " << i << endl;
      return false;
    }
  }
  return true;
}

bool Tester::testUpdateLocation() {
  DnaDb db(101, simpleHash, LINEAR);
  DNA d("GGTTA", 100000);//create dna obj                                                                                                                                         
  db.insert(d);//insert it                                                                                                                                                        
  if (!db.updateLocId(d, 100123)) return false;//update it loc                                                                                                                    
  if (db.getDNA("GGTTA", 100123).getLocId() != 100123) return false;//check if its updated                                                                                        
  return true;
}

bool Tester::testPolicySwitch() {
  DnaDb db(101, simpleHash, LINEAR);

  for (int i = 0; i < 40; ++i) {//insert dna with linear                                                                                                                          
    if (!db.insert(DNA("TTAA" + to_string(i), 100000 + i))) return false;
  }

  db.changeProbPolicy(QUADRATIC);//switch to quadratic                                                                                                                            

  for (int i = 40; i < 80; ++i) {//insert more                                                                                                                                    
    if (!db.insert(DNA("TTAA" + to_string(i), 100000 + i))) return false;
  }

  for (int i = 0; i < 30; ++i) db.rehashIfNeeded();//rehash                                                                                                                       

  for (int i = 0; i < 80; ++i) {//check all are still accessible                                                                                                                  
    DNA result = db.getDNA("TTAA" + to_string(i), 100000 + i);
    if (result.getSequence().empty()) {
      cout << "GetDNA failed for TTAA" << i << endl;
      return false;
    }
  }
  return true;
}

bool Tester::testFindFromOldTable() {
  DnaDb db(101, simpleHash, LINEAR);
  for (int i = 0; i < 60; ++i) {//insert 60 dna to pass 0.5                                                                                                                       
    if (!db.insert(DNA("CCGG" + to_string(i), 100000 + i))) return false;
  }

  for (int i = 0; i < 30; ++i) db.rehashIfNeeded();//rehash to move some                                                                                                          
  //make sure they are all still accessible                                                                                                                                       
  for (int i = 0; i < 60; ++i) {
    DNA result = db.getDNA("CCGG" + to_string(i), 100000 + i);
    if (result.getSequence().empty()) {
      cout << "GetDNA failed for CCGG" << i << endl;
      return false;
    }
  }
  return true;
}

int main() {
  Tester tester;
  cout << "testInsertNoCollision: " << endl;
  if(tester.testInsertNoCollision()) cout << "pass" << endl;
  else cout<< "fail" << endl;
  cout << "testInsertWithCollision: " << endl;
  if(tester.testInsertWithCollision()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testRejectDuplicate: " << endl;
  if(tester.testRejectDuplicate()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testInvalidLocation: " << endl;
  if(tester.testInvalidLocation()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testRemoveSuccess: " << endl;
  if(tester.testRemoveSuccess()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testRemoveNonExistent: " << endl;
  if(tester.testRemoveNonExistent()) cout<< "pass" << endl;
  else cout << "fail" << endl;
  cout << "testRehashByLoadFactor: " << endl;
  if(tester.testRehashByLoadFactor()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testRehashByDeleteRatio: " << endl;
  if(tester.testRehashByDeleteRatio()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testUpdateLocation: " << endl;
  if(tester.testUpdateLocation()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testPolicySwitch: " << endl;
  if(tester.testPolicySwitch()) cout << "pass" << endl;
  else cout << "fail" << endl;
  cout << "testFindFromOldTable: " << endl;
  if(tester.testFindFromOldTable()) cout << "pass" << endl;
  else cout << "fail" << endl;
  return 0;
}



