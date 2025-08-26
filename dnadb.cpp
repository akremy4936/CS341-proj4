File Edit Options Buffers Tools C++ Help                                                                                                                                          
// CMSC 341 - Spring 25 - Project 4                                                                                                                                               
#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    m_hash = hash;
    m_currProbing = probing;
    m_newPolicy = probing;
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    m_transferIndex = 0;

    //Set current capacity with prime bounds                                                                                                                                      
    if (size < MINPRIME)
        m_currentCap = MINPRIME;
    else if (size > MAXPRIME)
        m_currentCap = MAXPRIME;
    else if (!isPrime(size))
        m_currentCap = findNextPrime(size);
    else
        m_currentCap = size;

    //Allocate and initialize table                                                                                                                                               
    m_currentTable = new DNA*[m_currentCap];
    for (int i = 0; i < m_currentCap; i++)
        m_currentTable[i] = nullptr;

    m_currentSize = 0;
    m_currNumDeleted = 0;
}

DnaDb::~DnaDb(){//destructor                                                                                                                                                      
    for (int i = 0; i < m_currentCap; i++)
        delete m_currentTable[i];
    delete[] m_currentTable;

    if (m_oldTable != nullptr) {
        for (int i = 0; i < m_oldCap; i++)
            delete m_oldTable[i];
        delete[] m_oldTable;
    }
}

void DnaDb::changeProbPolicy(prob_t policy){//changes the probing policy(linear, quadratic, double hash)                                                                          
    m_newPolicy = policy;
}
bool DnaDb::insert(DNA dna) {//insert dna into current table                                                                                                                      
  if (dna.getLocId() < MINLOCID || dna.getLocId() > MAXLOCID)//wont insert if not correct location                                                                                
    return false;
  if (!getDNA(dna.getSequence(), dna.getLocId()).getSequence().empty())//wont take duplicates                                                                                     
    return false;

  int key = m_hash(dna.getSequence());
  int index = key % m_currentCap;
  int i = 0;

  while (i < m_currentCap) {//probe index for deleted/empty spot                                                                                                                  
    int probe;
    if (m_currProbing == LINEAR)
      probe = (index + i) % m_currentCap;
    else if (m_currProbing == QUADRATIC)
      probe = (index + i * i) % m_currentCap;
    else
      probe = (index + i * (11 - (key % 11))) % m_currentCap;

    if (m_currentTable[probe] == nullptr) {//if not used spot, insert                                                                                                             
      m_currentTable[probe] = new DNA(dna);
      m_currentTable[probe]->setUsed(true);
      m_currentSize++;
      rehashIfNeeded();
      return true;
    }
    else if (!m_currentTable[probe]->getUsed()) {//if been deleted, use for new data                                                                                              
      *m_currentTable[probe] = dna;
      m_currentTable[probe]->setUsed(true);
      m_currNumDeleted--;
      rehashIfNeeded();
      return true;
    }
    i++;
  }
  return false;
}

bool DnaDb::remove(DNA dna) {
    int key = m_hash(dna.getSequence());
    int index = key % m_currentCap;
    int i = 0;

    //Try to remove from current table                                                                                                                                            
    while (i < m_currentCap) {
        int probe;
        if (m_currProbing == LINEAR)
            probe = (index + i) % m_currentCap;
        else if (m_currProbing == QUADRATIC)
            probe = (index + i * i) % m_currentCap;
        else
            probe = (index + i * (11 - (key % 11))) % m_currentCap;

        if (m_currentTable[probe] == nullptr)
            break;

        if (*m_currentTable[probe] == dna && m_currentTable[probe]->getUsed()) {
            m_currentTable[probe]->setUsed(false);
            m_currNumDeleted++;

            //Transfer 25% if rehash is in progress OR start rehash if needed                                                                                                     
            rehashIfNeeded();
            return true;
        }

        i++;
    }

    //Try to remove from old table if it exists                                                                                                                                   
    if (m_oldTable != nullptr) {
        key = m_hash(dna.getSequence());
        index = key % m_oldCap;
        i = 0;

        while (i < m_oldCap) {
            int probe;
            if (m_oldProbing == LINEAR)
                probe = (index + i) % m_oldCap;
            else if (m_oldProbing == QUADRATIC)
                probe = (index + i * i) % m_oldCap;
            else
                probe = (index + i * (11 - (key % 11))) % m_oldCap;

            if (m_oldTable[probe] == nullptr)
                break;

            if (*m_oldTable[probe] == dna && m_oldTable[probe]->getUsed()) {
                m_oldTable[probe]->setUsed(false);
                m_oldNumDeleted++;

                //Continue incremental transfer if already rehashing                                                                                                              
                rehashIfNeeded();
                return true;
            }

            i++;
        }
    }

    return false; // Not found                                                                                                                                                    
}

const DNA DnaDb::getDNA(string sequence, int location) const {
    DNA target(sequence, location, true);

    //Search through the current table                                                                                                                                            
    int key = m_hash(sequence);
    int index = key % m_currentCap;
    int i = 0;

    while (i < m_currentCap) {
        int probe;
        if (m_currProbing == LINEAR)
            probe = (index + i) % m_currentCap;
        else if (m_currProbing == QUADRATIC)
            probe = (index + i * i) % m_currentCap;
        else
            probe = (index + i * (11 - (key % 11))) % m_currentCap;

        //Only break if the spot is empty and it was never used                                                                                                                   
        if (m_currentTable[probe] == nullptr)
            break;

        if (*m_currentTable[probe] == target)
            return *m_currentTable[probe];

        i++;
    }

    //Search through old table                                                                                                                                                    
    if (m_oldTable != nullptr) {
        key = m_hash(sequence);
        index = key % m_oldCap;
        i = 0;

        while (i < m_oldCap) {
            int probe;
            if (m_oldProbing == LINEAR)
                probe = (index + i) % m_oldCap;
            else if (m_oldProbing == QUADRATIC)
                probe = (index + i * i) % m_oldCap;
            else
                probe = (index + i * (11 - (key % 11))) % m_oldCap;

            if (m_oldTable[probe] == nullptr)
                break;

            if (*m_oldTable[probe] == target)
                return *m_oldTable[probe];

            i++;
        }
    }

    return DNA();
}

bool DnaDb::updateLocId(DNA dna, int location){
  DNA result = getDNA(dna.getSequence(), dna.getLocId());//check if it exists                                                                                                     
  if (result.getSequence().empty())
    return false;

  //Update in current table                                                                                                                                                       
  int key = m_hash(dna.getSequence());
  int index = key % m_currentCap;
  int i = 0;

  while (i < m_currentCap) {
    int probe;
    if (m_currProbing == LINEAR)
      probe = (index + i) % m_currentCap;
    else if (m_currProbing == QUADRATIC)
      probe = (index + i * i) % m_currentCap;
    else
      probe = (index + i * (11 - (key % 11))) % m_currentCap;
    //object not here if bucket empty                                                                                                                                             
    if (m_currentTable[probe] == nullptr)
      break;
    if (*m_currentTable[probe] == dna) {
      m_currentTable[probe]->setLocID(location);
      return true;
    }
    i++;
  }

  //Update in old table                                                                                                                                                           
  if (m_oldTable != nullptr) {
    key = m_hash(dna.getSequence());
    index = key % m_oldCap;
    i = 0;

    while (i < m_oldCap) {
      int probe;
      if (m_oldProbing == LINEAR)
        probe = (index + i) % m_oldCap;
      else if (m_oldProbing == QUADRATIC)
        probe = (index + i * i) % m_oldCap;
      else
        probe = (index + i * (11 - (key % 11))) % m_oldCap;

      if (m_oldTable[probe] == nullptr)
        break;

      if (*m_oldTable[probe] == dna) {
        m_oldTable[probe]->setLocID(location);
        return true;
      }

      i++;
    }
  }
  return false;

}

float DnaDb::lambda() const {
    return (float)m_currentSize / m_currentCap;
}

float DnaDb::deletedRatio() const {
    if (m_currentSize == 0) return 0.0;
    return (float)m_currNumDeleted / m_currentSize;
}

void DnaDb::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}
bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}
int DnaDb::findNextPrime(int current){
    //always stay within the range [MINPRIME-MAXPRIME]                                                                                                                            
    //the smallest prime starts at MINPRIME                                                                                                                                       
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) {
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0)
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME                                                                                                                                         
    return MAXPRIME;
}
void DnaDb::rehashIfNeeded() {
  //live inserts                                                                                                                                                                  
  int liveCount = m_currentSize - m_currNumDeleted;

  // Check whether we need to start rehashing                                                                                                                                     
  if (m_oldTable == nullptr && (lambda() > 0.5 || deletedRatio() > 0.8)) {
    //Save old table info                                                                                                                                                         
    m_oldTable = m_currentTable;
    m_oldCap = m_currentCap;
    m_oldSize = m_currentSize;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldProbing = m_currProbing;

    //Compute new capacity                                                                                                                                                        
    m_currentCap = findNextPrime(liveCount * 4);
    m_currentTable = new DNA*[m_currentCap];
    for (int i = 0; i < m_currentCap; i++) {
      m_currentTable[i] = nullptr;
    }
    // Reset new table stats                                                                                                                                                      
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_currProbing = m_newPolicy;

    m_transferIndex = 0;
    //Begin transferring                                                                                                                                                          
    transferData();
  }
  //If already rehashing, continue transfer                                                                                                                                       
  else if (m_oldTable != nullptr) {
    transferData();
  }
}

bool DnaDb::insertNoRehash(const DNA& dna) {
  int key = m_hash(dna.getSequence());
  int index = key % m_currentCap;
  int i = 0;

  while (i < m_currentCap) {
    int probe;
    if (m_currProbing == LINEAR)
      probe = (index + i) % m_currentCap;
    else if (m_currProbing == QUADRATIC)
      probe = (index + i * i) % m_currentCap;
    else
      probe = (index + i * (11 - (key % 11))) % m_currentCap;

    if (m_currentTable[probe] == nullptr) {//insert into first found nullptr bucket                                                                                               
      m_currentTable[probe] = new DNA(dna);
      m_currentTable[probe]->setUsed(true);
      m_currentSize++;
      return true;
    }
    i++;
  }
  return false; //should not happen if rehash size is correct                                                                                                                     
}

void DnaDb::transferData() {
  if (m_oldTable == nullptr) return;
  //find number that has be transferred                                                                                                                                           
  int chunkSize = m_oldCap / 4;
  int end = m_transferIndex + chunkSize;
  if (end > m_oldCap) end = m_oldCap;
  //go through old table                                                                                                                                                          
  for (int i = m_transferIndex; i < end; ++i) {
    if (m_oldTable[i] != nullptr && m_oldTable[i]->getUsed()) {//insert non deleted into the new table                                                                            
      insertNoRehash(*m_oldTable[i]);  //avoids triggering rehash again                                                                                                           
    }

    delete m_oldTable[i];
    m_oldTable[i] = nullptr;
  }

  m_transferIndex = end;
  //when done transferring clean the old table                                                                                                                                    
  if (m_transferIndex >= m_oldCap) {
    delete[] m_oldTable;
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
  }
}



