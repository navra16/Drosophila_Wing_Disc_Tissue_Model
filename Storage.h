#ifndef STORAGE_H_
#define STORAGE_H_

#include <fstream>	     
#include <memory>
#include <iomanip>

class System;

class Storage {
	
	
	std::weak_ptr<System> system;
	//std::shared_ptr<ExternalForce> grip;
	std::ofstream output;
	std::ofstream statesOutput;
	
	
	int iteration = 0;
	int iteration_daughter = 0;
	int iteration2 = 0;

public: 
	Storage(std::weak_ptr<System>);

	void storeVariables(void);

	void print_VTK_File(void);


};

#endif /*STORAGE_H_*/