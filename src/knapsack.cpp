#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <algorithm>

typedef std::vector<std::pair<int, int> > ItemVector;

struct inst
{
	ItemVector items; //weight, cost
	int id;
	int max_weight;
	int num_items;
};
struct soln
{
	uint64_t config;
	int id;
	int num_items;
	int cost;
};

int bbTopCost;
std::vector<int> bbTopConfig;

//helpers
void printInst(inst instance)
{
	printf("Id: %d, Max Weight: %d, Num Items: %d\n", instance.id, instance.max_weight, instance.num_items);
	for (int i = 0; i < instance.num_items; i++)
	{
		printf("Item %d: Weight=%d, Cost=%d\n", i, instance.items[i].first, instance.items[i].second);
	}
	printf("\n");
}

void printSoln(soln solution)
{
	printf("Id: %d, n: %d, cost: %d, config: %lld\n", solution.id, solution.num_items, solution.cost, solution.config);
}

int GetInstances(std::vector<inst> &instances, const char* filename, bool verbose) 
{
	std::string line;
  	std::ifstream infile(filename);

  	if (infile)
  	{
	  	while (std::getline(infile, line))
	  	{
  			inst newInst;
	  		std::istringstream iss(line);

	  		iss >> newInst.id;
	  		iss >> newInst.num_items;
	  		iss >> newInst.max_weight;
	  		newInst.items.resize(newInst.num_items);

	  		for (int i = 0; i <= newInst.num_items; i++)
	  		{
	  			iss >> newInst.items[i].first;
	  			iss >> newInst.items[i].second;
	  		}

	  		instances.push_back(newInst);

	  		if (verbose)
	  		{
	  			printf("Pushed New Instance:\n");
	  			printInst(newInst);
	  		}
	  	}
  	}

  	return 0;
}

int DropCostBits(std::vector<inst> &instances, int numBits, bool verbose)
{
	if (numBits == 0) return 0;
	for (int i = 0; i < instances.size(); i++)
	{
		if (verbose) printf("Assigning instance %d items to bins...\n", i);

		for (int j = 0; j < instances[i].num_items; j++)
		{
			instances[i].items[j].second = instances[i].items[j].second & ~((int)std::pow(2, numBits) - 1);
		}
	}
	return 0;
}

int RoundCosts(std::vector<inst> &instances, int binSize, bool verbose)
{
	if (binSize == 1) return 0;
	for (int i = 0; i < instances.size(); i++)
	{
		if (verbose) printf("Assigning instance %d items to bins...\n", i);

		for (int j = 0; j < instances[i].num_items; j++)
		{
			int margin = binSize - (instances[i].items[j].second % binSize);
			instances[i].items[j].second += margin;
			if (verbose) printf("Rounding item %d cost up by %d\n", j, margin);
		}
	}
	return 0;
}

//algos
int Bruteforce(std::vector<inst> instances, std::vector<soln> &solutions, bool verbose)
{
	//start timer
	std::clock_t start;
    double duration;

    start = std::clock();

	while (!instances.empty())
	{
		inst curInst = instances.back();
		instances.pop_back();

		int topCost = 0;
		uint64_t topConfig = 0;

		for (uint64_t i = 0; i < std::pow(2, curInst.num_items); i++)
		{
			//if (verbose) printf("Attempting config: %lld\n", i);
			int curCost = 0;
			int curWeight = 0;
			std::bitset<64> curConfig(i);

			for (int j = 0; j <= curInst.num_items; j++)
			{
				if (curConfig[j])
				{
					curCost += curInst.items[j].second;
					curWeight += curInst.items[j].first;
					//if (verbose) printf("Adding item. Weight: %d, Cost: %d\n", curInst.items[j].first, curInst.items[j].second);
				}
			}
			if (curCost > topCost && curWeight <= curInst.max_weight)
			{
				//record potential solution
				topCost = curCost;
				topConfig = i;
			}
		}

		soln newSoln = {
			.config 	= topConfig,
			.id 		= curInst.id,
			.num_items 	= curInst.num_items,
			.cost 		= topCost
		};
		solutions.push_back(newSoln);
		if (verbose) printSoln(newSoln);
	}

	//stop timer
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("%f\n", duration);

	return 0;
}

int Greeeeeedy(std::vector<inst> instances, std::vector<soln> &solutions, bool verbose)
{
	//start timer
	std::clock_t start;
    double duration;

    start = std::clock();

    //process
	while (!instances.empty())
	{
		inst curInst = instances.back();
		instances.pop_back();

		int solnWeight = 0;
		int solnCost = 0;
		std::bitset<40> solnConfig;

		for (int j = 0; j < curInst.num_items; j++)
		{
			int topCWR = 0;
			int topCWRindex = -1;

			for (int i = 0; i < curInst.num_items; i++)//find highest CWR
			{
				int curCWR = curInst.items[i].second / curInst.items[i].first;
				if (curCWR > topCWR && solnWeight + curInst.items[i].first <= curInst.max_weight && !solnConfig.test(i))
				{
					topCWR = curCWR;
					topCWRindex = i;
				}
			}
			if (topCWRindex != -1)//add it to the solution
			{
				solnWeight += curInst.items[topCWRindex].first;
				solnCost += curInst.items[topCWRindex].second;
				solnConfig.set(topCWRindex);
			}
		}

		soln newSoln = {
			.config 	= solnConfig.to_ullong(),
			.id 		= curInst.id,
			.num_items 	= curInst.num_items,
			.cost 		= solnCost
		};
		solutions.push_back(newSoln);
		if (verbose) printSoln(newSoln);
	}

	//stop timer
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("%f\n", duration);

	return 0;
}

int BBRecursion(ItemVector items, std::vector<int> config, int max_weight)
{
	int firstUndecidedIndex = -1;
	int includedCost = 0;
	int includedWeight = 0;
	int undecidedCost = 0;
	for (int i = 0; i < config.size(); ++i)
	{
		if (config[i] == 1)
		{
			includedCost += items[i].second;
			includedWeight += items[i].first;
		}
		else if (config[i] == -1)
		{
			undecidedCost += items[i].second;

			if (firstUndecidedIndex == -1)
			{
				firstUndecidedIndex = i;
			}
		}
	}

	if (includedWeight <= max_weight)
	{
		if (includedCost > bbTopCost)
		{
			bbTopCost = includedCost;
			bbTopConfig = config;
		}

		if (includedCost + undecidedCost > bbTopCost) //potential better solution exists
		{
			//recursive call using undecidedIndex
			std::vector<int> includeNext(config);
			includeNext[firstUndecidedIndex] = 1;
			std::vector<int> excludeNext(config);
			excludeNext[firstUndecidedIndex] = 0;
			BBRecursion(items, includeNext, max_weight);
			BBRecursion(items, excludeNext, max_weight);
		}
	}
	return 0;
}

int BBWrapper(std::vector<inst> instances, std::vector<soln> &solutions, bool verbose)
{
	std::clock_t start;
    double duration;

    start = std::clock();

    //process
	while (!instances.empty())
	{
		//printf("New Instance\n");
		inst curInst = instances.back();
		instances.pop_back();

		bbTopConfig.clear();
		bbTopCost = 0;

		//run recursion
		std::vector<int> config (curInst.num_items, -1);//initialize config vector// -1: undecided, 0: excluded, 1: included
		BBRecursion(curInst.items, config, curInst.max_weight);

		std::bitset<64> solnConfig;
		if (verbose) {
			for (int i = 0; i < curInst.num_items; ++i)
			{
				if (bbTopConfig[i] == 1)
				{
					solnConfig[i] = 1;
				}
			}
		}

		soln newSoln = {
			.config 	= solnConfig.to_ullong(),
			.id 		= curInst.id,
			.num_items 	= curInst.num_items,
			.cost 		= bbTopCost
		};
		solutions.push_back(newSoln);
		if (verbose) printSoln(newSoln);
	}

	//stop timer
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("%f\n", duration);

	return 0;
}

int Dynamic(std::vector<inst> instances, std::vector<soln> &solutions, bool verbose)
{
	std::clock_t start;
    double duration;

    start = std::clock();

	while (!instances.empty())
	{
		inst curInst = instances.back();
		instances.pop_back();

		std::vector<std::vector<int> > costs(curInst.num_items + 1, std::vector<int> (curInst.max_weight + 1, 0));
		for (int i = 1; i <= curInst.num_items; i++)
		{
			for (int j = 0; j <= curInst.max_weight; j++)
			{
				if (j > curInst.items[i - 1].first)
				{
					costs[i][j] = std::max(
						costs[i - 1][j], 
						costs[i - 1][j - curInst.items[i - 1].first] + curInst.items[i - 1].second
					);
				}
				else
				{
					costs[i][j] = costs[i - 1][j];
				}
			}
		}

		soln newSoln = {
			.config 	= 0,
			.id 		= curInst.id,
			.num_items 	= curInst.num_items,
			.cost 		= costs[curInst.num_items][curInst.max_weight]
		};
		solutions.push_back(newSoln);
		if (verbose) printSoln(newSoln);
	}

	//stop timer
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("%f\n", duration);

	return 0;
}

int FPTAS(std::vector<inst> instances, std::vector<soln> &solutions, int costStepSize, bool verbose)
{
	std::clock_t start;
    double duration;

    start = std::clock();

	while (!instances.empty())
	{
		inst curInst = instances.back();
		instances.pop_back();

		std::vector<std::vector<int> > costs(curInst.num_items + 1, std::vector<int> (curInst.max_weight + 1, 0));
		for (int i = 1; i <= curInst.num_items / costStepSize; i++)
		{
			for (int j = 0; j <= curInst.max_weight; j++)
			{
				if (j > curInst.items[i - 1].first)
				{
					costs[i][j] = std::max(
						costs[i - 1][j], 
						costs[i - 1][j - curInst.items[i - 1].first] + curInst.items[i - 1].second
					);
				}
				else
				{
					costs[i][j] = costs[i - 1][j];
				}
			}
		}

		soln newSoln = {
			.config 	= 0,
			.id 		= curInst.id,
			.num_items 	= curInst.num_items,
			.cost 		= costs[curInst.num_items / costStepSize][curInst.max_weight]
		};
		solutions.push_back(newSoln);
		if (verbose) printSoln(newSoln);
	}

	//stop timer
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("%f\n", duration);

	return 0;
}

int main(int argc, const char** argv) 
{
	if (argc < 3 || argc > 5){
		printf("Error: Unexpected argument count.\nUsage: <filename> <algorithm> [FPTAS num ignored bits] [-v]\n\n");
		printf("ALGORITHM:\n1 - bruteforce\n2 - Greedy (cost/weight ratio heuristic)\n");
		printf("3 - Branch and Bound (cost as bounding factor)\n4 - Dynamic Programming (cost decomposition)\n");
		printf("5 - FPTAS\n");
		return 1;
	}

	std::vector<inst> instances;
	std::vector<soln> solutions;

	bool verbose;
	if (argc == 4)
	{
		verbose = (std::strncmp(argv[3], "-v", 2) == 0) ? 1 : 0;
	}
	else if(argc == 5)
	{
		verbose = (std::strncmp(argv[4], "-v", 2) == 0) ? 1 : 0;
	}

	GetInstances(instances, argv[1], verbose);

	int algo;

	try //C++ style!
	{
		algo = std::stoi(argv[2]);
	}
	catch(...)
	{
		printf("Error: Invalid algorithm selection");
		return 1;
	}

	switch (algo)
	{
		case 1:
			printf("Running bruteforce algorithm...\n");
			Bruteforce(instances, solutions, verbose);
			break;
		case 2:
			printf("Running greedy algorithm...\n");
			Greeeeeedy(instances, solutions, verbose);
			break;
		case 3:
			printf("Running recursive B&B algorithm...\n");
			BBWrapper(instances, solutions, verbose);
			break;
		case 4:
			printf("Running dynamic programming algorithm...\n");
			Dynamic(instances, solutions, verbose);
			break;
		case 5:
			if (argc < 4)
			{
				printf("Error: Invalid bin size.\n");
				return 1;
			}
			int errorLevel = atoi(argv[3]); // C style!
			if (!errorLevel)
			{
				printf("Error: Invalid bin size.\n");
				return 1;
			}
			printf("Running FPTAS algorithm (num cost bits ignored: %d)...\n", errorLevel);
			RoundCosts(instances, errorLevel, verbose);
			FPTAS(instances, solutions, errorLevel, verbose);
			break;
	}

	instances.clear();
	solutions.clear();
	bbTopConfig.clear();

	return 0;
}



