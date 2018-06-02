#include "Graph.h"
#include <chrono>
#include <iostream>

using namespace std;

int main()
{
	Graph graph = Graph();
	graph.readGraph("PrimaTest.txt");
	system("pause");
	return 0;
}