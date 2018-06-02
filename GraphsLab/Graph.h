#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <queue>
#include <functional>
#include <stack>
#include <iostream>

using namespace std;

class AdjMatrixGraphRepresentation;
class AdjListGraphRepresentation;
class EdgeListGraphRepresentation;
class Graph;
class Edge;

typedef vector<vector<Edge*>*>* AdjacencyMatrix;
typedef vector<map<int, Edge*>*>* AdjacencyList;
typedef map<pair<int, int>, Edge*>* EdgeList;
typedef map<int, int>* Vertexes;
typedef vector<map<int, bool>*>* MarkedEdges;


class StringParser
{
public:
	static vector<string> split(const string& str, char c)
	{
		vector<string> result;
		istringstream stream(str);
		string token;
		while (getline(stream, token, c))
			result.push_back(token);
		return result;
	}
};

class DSU {
private:
	vector<int>* parent;
	vector<int>* size;

public:
	DSU() {}
	DSU(int n) {
		parent = new vector<int>(n);
		size = new vector<int>(n);
	}
	~DSU() {
		delete parent;
		delete size;
	}

	void makeSet(int v) {
		(*parent)[v] = v;
		(*size)[v] = 1;
	}

	int find(int x) {
		if ((*parent)[x] == x)
			return x;
		int xRoot = find((*parent)[x]);
		(*parent)[x] = xRoot;
		return xRoot;
	}

	void unite(int x, int y) {
		int xRoot = find(x);
		int yRoot = find(y);
		if (xRoot != yRoot) {
			if ((*size)[xRoot] >= (*size)[yRoot])
			{
				(*size)[xRoot] += (*size)[yRoot];
				(*parent)[yRoot] = xRoot;
			}
			else
			{
				(*size)[yRoot] += (*size)[xRoot];
				(*parent)[xRoot] = yRoot;
			}
		}
	}

	map<int, int> getRoots() {
		map<int, int> roots = map<int, int>();
		for (int i = 0; i < parent->size(); ++i)
			if (i == (*parent)[i])
				roots[i] = (*size)[i];
		return roots;
	}
};


class Edge {
	int pv, pw, pcap, pflow;
	bool duplicate;
public:
	Edge(int v, int w, int cap) : pv(v), pw(w), pcap(cap), pflow(0), duplicate(false) {}
	Edge(int v, int w, int cap, bool isDuplicate) : pv(v), pw(w), pcap(cap), pflow(0), duplicate(isDuplicate) {}
	int v() const { return pv; }
	int w() const { return pw; }
	int cap() const { return pcap; }
	int flow() const { return pflow; }
	bool from(int v) const { return pv == v; }
	int other(int v) const {
		return from(v) ? pw : pv;
	}
	int capRto(int v, int d) {
		return from(v) ? pflow : pcap - pflow;
	}
	void addflowRto(int v, int d) {
		pflow += from(v) ? -d : d;
	}
	void setCap(int newCap) {
		pcap = newCap;
	}
	bool isDuplicate() {
		return duplicate;
	}
	int flowLeft() {
		return cap() - flow();
	}
};

class GraphRepresentation
{
protected:
	bool isOriented, isWeighted;
	Vertexes vertexDegrees = nullptr;
	virtual DSU* buildDSU() = 0;

public:

	int edgesAmount = 0;
	GraphRepresentation() = default;
	virtual ~GraphRepresentation() = default;

	virtual void readGraph(istream& stream, vector<string> parameters) = 0;
	virtual void addEdge(int from, int to, int weight) = 0;
	virtual void removeEdge(int from, int to) = 0;
	virtual int changeEdge(int from, int to, int newWeight) = 0;
	virtual AdjMatrixGraphRepresentation* transformToAdjMatrix() = 0;
	virtual AdjListGraphRepresentation* transformToAdjList() = 0;
	virtual EdgeListGraphRepresentation* transformToListOfEdges() = 0;
	virtual void writeGraph(string fileName) = 0;
	virtual bool isAdjMatrixGraph() = 0;
	virtual bool isAdjListGraph() = 0;
	virtual bool isEdgeListGraph() = 0;
	int getBeginVertex(bool& circleExists) {
		int evenVertexesAmount = 0;
		int oddVertex = -1;
		int evenVertex = -1;
		for (auto it = vertexDegrees->begin(); it != vertexDegrees->end(); ++it) {
			if (it->second % 2 == 0)
			{
				++evenVertexesAmount;
				if (evenVertex == -1 && it->second > 0)
					evenVertex = it->first;
			}
			else
				oddVertex = it->first;
		}
		circleExists = evenVertexesAmount == vertexDegrees->size();
		if (circleExists)
			return evenVertex;
		if (vertexDegrees->size() - evenVertexesAmount <= 2)
			return oddVertex;
		return -1;
	}
	bool connectedComponentsCheck() {
		DSU* dsu = buildDSU();
		map<int, int> roots = dsu->getRoots();
		if (roots.size() <= 1)
			return true;
		int nonEmptyComponentsAmount = 0;
		for (auto it = roots.begin(); it != roots.end(); ++it)
			if (it->second > 1)
				nonEmptyComponentsAmount++;
		return nonEmptyComponentsAmount <= 1;
	}

	virtual int calcWeight() = 0;
	virtual int calcFlow() = 0;

	friend Graph;
};

class AdjMatrixGraphRepresentation : GraphRepresentation
{
	AdjacencyMatrix adjacencyMatrix = nullptr;
	AdjacencyMatrix readMatrix(istream& stream, int vertexAmount)
	{
		const AdjacencyMatrix matrix = new vector<vector<Edge*>*>(vertexAmount);
		vertexDegrees = new map<int, int>();
		for (size_t i = 0; i < vertexAmount; ++i)
		{
			(*vertexDegrees)[i] = 0;
			(*matrix)[i] = new vector<Edge*>(vertexAmount);
		}
		for (size_t i = 0; i < vertexAmount; ++i)
		{
			string line;
			getline(stream, line);
			vector<string> tokens = StringParser::split(line, ' ');
			for (size_t j = 0; j < vertexAmount; ++j) {
				const int weight = stoi(tokens[j]);
				if (weight > 0 && matrix->at(i)->at(j) == nullptr && (matrix->at(j)->at(i) == nullptr || isOriented)) {
					(*(*matrix)[i])[j] = new Edge(i, j, weight);
					++(*vertexDegrees)[i];
					++edgesAmount;
					if (!isOriented)
						++(*vertexDegrees)[j];
				}
			}
		}
		return matrix;
	}
	static void writeMatrix(ostream& stream, AdjacencyMatrix matrix, bool isOriented)
	{
		for (size_t i = 0; i < matrix->size(); ++i)
		{
			for (size_t j = 0; j < matrix->size(); ++j)
			{
				if (matrix->at(i)->at(j) != nullptr)
					stream << (*(*matrix)[i])[j]->cap();
				else
					stream << "0";
				if (j < matrix->size() - 1)
					stream << " ";
			}
			if (i < matrix->size() - 1)
				stream << endl;
		}
	}
	void copy(const AdjMatrixGraphRepresentation& graph) {
		this->isOriented = graph.isOriented;
		this->isWeighted = graph.isWeighted;
		this->edgesAmount = graph.edgesAmount;
		this->vertexDegrees = new map<int, int>(*graph.vertexDegrees);
		this->adjacencyMatrix = new vector<vector<Edge*>*>(graph.adjacencyMatrix->size());
		for (size_t i = 0; i < adjacencyMatrix->size(); ++i)
			(*adjacencyMatrix)[i] = new vector<Edge*>(*(*graph.adjacencyMatrix)[i]);
	}
	DSU* buildDSU() {
		DSU* dsu = new DSU(adjacencyMatrix->size());
		for (int i = 0; i < adjacencyMatrix->size(); ++i)
			dsu->makeSet(i);
		for (int i = 0; i < adjacencyMatrix->size(); ++i)
			for (int j = 0; j < adjacencyMatrix->size(); ++j)
				dsu->unite(i, j);
		return dsu;
	}

public:
	AdjMatrixGraphRepresentation() = default;
	AdjMatrixGraphRepresentation(int vertexAmount, int edgesAmount, bool isOriented, bool isWeighted)
	{
		this->vertexDegrees = new map<int, int>();
		this->adjacencyMatrix = new vector<vector<Edge*>*>(vertexAmount);
		for (int i = 0; i < vertexAmount; ++i)
			(*adjacencyMatrix)[i] = new vector<Edge*>(vertexAmount, 0);
		this->isOriented = isOriented;
		this->isWeighted = isWeighted;
		this->edgesAmount = edgesAmount;
	}
	AdjMatrixGraphRepresentation(const AdjMatrixGraphRepresentation& graph)
	{
		copy(graph);
	}
	~AdjMatrixGraphRepresentation()
	{
		delete vertexDegrees;
		for (int i = 0; i < adjacencyMatrix->size(); ++i) {
			for (int j = 0; j < adjacencyMatrix->at(i)->size(); ++j)
				delete(adjacencyMatrix->at(i)->at(j));
			delete (*adjacencyMatrix)[i];
		}
		delete adjacencyMatrix;
	};
	AdjMatrixGraphRepresentation& operator=(const AdjMatrixGraphRepresentation& graph) {
		copy(graph);
		return *this;
	}

	bool isAdjMatrixGraph() override { return true; }
	bool isAdjListGraph() override { return false; }
	bool isEdgeListGraph() override { return false; }

	void readGraph(istream& stream, vector<string> parameters) override
	{
		const int vertexAmount = stoi(parameters[1]);

		isOriented = parameters[2] == "1";
		isWeighted = parameters[3] == "1";

		adjacencyMatrix = readMatrix(stream, vertexAmount);
	}


	void addEdge(int from, int to, int weight) override
	{
		if (adjacencyMatrix->at(from)->at(to) == nullptr && (adjacencyMatrix->at(from)->at(to) == nullptr || isOriented)) {
			++edgesAmount;
			++(*vertexDegrees)[from];
			(*(*adjacencyMatrix)[from])[to] = new Edge(from, to, weight);
			if (!isOriented) {
				++(*vertexDegrees)[to];
				(*(*adjacencyMatrix)[to])[from] = new Edge(from, to, weight);
			}
		}
	}


	void removeEdge(int from, int to)  override
	{
		--edgesAmount;
		--(*vertexDegrees)[from];
		delete((*(*adjacencyMatrix)[from])[to]);
		(*(*adjacencyMatrix)[from])[to] = nullptr;
		if (!isOriented) {
			--(*vertexDegrees)[to];
			delete((*(*adjacencyMatrix)[to])[from]);
			(*(*adjacencyMatrix)[to])[from] = nullptr;
		}
	}


	int changeEdge(int from, int to, int newWeight)  override
	{
		const int oldWeight = (*(*adjacencyMatrix)[from])[to]->cap();
		(*(*adjacencyMatrix)[from])[to]->setCap(newWeight);
		return oldWeight;
	}


	AdjMatrixGraphRepresentation* transformToAdjMatrix()  override
	{
		return this;
	}


	AdjListGraphRepresentation* transformToAdjList()  override;


	EdgeListGraphRepresentation* transformToListOfEdges()  override;


	void writeGraph(string fileName) override
	{
		ofstream file(fileName);
		file << "C " << adjacencyMatrix->size() << endl;
		file << (isOriented ? "1 " : "0 ") << (isWeighted ? "1" : "0") << endl;

		writeMatrix(file, adjacencyMatrix, isOriented);

		file.close();
	}


	int calcWeight() override
	{
		return 0;
	}
	int calcFlow() override
	{
		return 0;
	}

	friend Graph;
	friend AdjListGraphRepresentation;
};

class AdjListGraphRepresentation : GraphRepresentation
{
	AdjacencyList adjacencyList = nullptr;
	AdjacencyList readList(istream& stream, int vertexAmount, bool isWeighted)
	{
		const AdjacencyList list = new vector<map<int, Edge*>*>(vertexAmount);
		vertexDegrees = new map<int, int>();
		for (size_t i = 0; i < vertexAmount; ++i)
		{
			(*vertexDegrees)[i] = 0;
			(*list)[i] = new map<int, Edge*>();
		}
		for (size_t i = 0; i < vertexAmount; ++i)
		{
			string line;
			getline(stream, line);
			vector<string> tokens = StringParser::split(line, ' ');
			for (size_t j = 0; j < tokens.size(); ++j)
			{
				int vertexId, weight = 1;
				if (isWeighted)
				{
					vertexId = stoi(tokens[j]) - 1;
					weight = stoi(tokens[++j]);
				}
				else
					vertexId = stoi(tokens[j]) - 1;
				if (list->at(i)->count(vertexId) == 0 && (list->at(vertexId)->count(i) == 0 || isOriented))
				{
					++edgesAmount;
					++(*vertexDegrees)[i];
					(*(*list)[i])[vertexId] = new Edge(i, vertexId, weight);
					if (!isOriented)
						++(*vertexDegrees)[vertexId];
				}
			}
		}
		return list;
	}
	static void writeList(ostream& stream, AdjacencyList list, bool isWeighted)
	{
		for (size_t i = 0; i < list->size(); ++i)
		{
			map<int, Edge*>* currentVertexMap = (*list)[i];
			for (auto it = currentVertexMap->begin(); it != currentVertexMap->end(); ++it)
			{
				stream << it->second->w() + 1;
				if (isWeighted)
					stream << " " << it->second->cap();
				if (next(it) != currentVertexMap->end())
					stream << " ";
			}
			if (i < list->size() - 1)
				stream << endl;
		}
	}
	void removeEdgeSimplex(AdjacencyList adjacencyList, int from, int to)
	{
		auto it = adjacencyList->at(from)->begin();
		bool edgeFound = false;
		for (; it != adjacencyList->at(from)->end(); ++it) {
			if (it->second->w() == to) {
				edgeFound = true;
				break;
			}
		}
		if (edgeFound) {
			delete(it->second);
			(*(*adjacencyList)[from]).erase(it);
			--(*vertexDegrees)[from];
		}
	}
	void copy(const AdjListGraphRepresentation& graph) {
		this->isOriented = graph.isOriented;
		this->isWeighted = graph.isWeighted;
		this->edgesAmount = graph.edgesAmount;
		this->vertexDegrees = new map<int, int>(*graph.vertexDegrees);
		this->adjacencyList = new vector<map<int, Edge*>*>(graph.adjacencyList->size());
		for (size_t i = 0; i < adjacencyList->size(); ++i)
			(*adjacencyList)[i] = new map<int, Edge*>(*(*graph.adjacencyList)[i]);
	}
	DSU* buildDSU() {
		DSU* dsu = new DSU(adjacencyList->size());
		for (int i = 0; i < adjacencyList->size(); ++i)
			dsu->makeSet(i);
		for (int i = 0; i < adjacencyList->size(); ++i) {
			map<int, Edge*>* edges = (*adjacencyList)[i];
			for (auto it = edges->begin(); it != edges->end(); ++it)
				dsu->unite(i, it->second->w());
		}
		return dsu;
	}

public:
	AdjListGraphRepresentation() = default;
	AdjListGraphRepresentation(int vertexAmount, int edgesAmount, bool isOriented, bool isWeighted)
	{
		this->vertexDegrees = new map<int, int>();
		this->edgesAmount = edgesAmount;
		this->adjacencyList = new vector<map<int, Edge*>*>(vertexAmount);
		for (int i = 0; i < adjacencyList->size(); ++i)
			(*adjacencyList)[i] = new map<int, Edge*>();
		this->isOriented = isOriented;
		this->isWeighted = isWeighted;
	}
	AdjListGraphRepresentation(int vertexAmount, int edgesAmount, bool isOriented, bool isWeighted, AdjacencyList adjList)
		: AdjListGraphRepresentation(vertexAmount, edgesAmount, isOriented, isWeighted) {
		for (size_t i = 0; i < adjList->size(); ++i)
			for (auto it = adjList->at(i)->begin(); it != adjList->at(i)->end(); ++it) {
				addEdge(i, it->first, it->second->cap());
				adjacencyList->at(i)->at(it->second->w())->addflowRto(it->second->w(), it->second->flow());
			}
	}
	AdjListGraphRepresentation(const AdjListGraphRepresentation& graph) {
		copy(graph);
	}
	~AdjListGraphRepresentation()
	{
		delete vertexDegrees;
		for (size_t i = 0; i < adjacencyList->size(); ++i) {
			for (auto it = adjacencyList->at(i)->begin(); it != adjacencyList->at(i)->end(); ++it)
				delete(it->second);
			delete (*adjacencyList)[i];
		}
		delete adjacencyList;
	}
	AdjListGraphRepresentation& operator=(const AdjListGraphRepresentation& graph) {
		copy(graph);
		return *this;
	}

	bool isAdjMatrixGraph() override { return false; }
	bool isAdjListGraph() override { return true; }
	bool isEdgeListGraph() override { return false; }

	void readGraph(istream& stream, vector<string> parameters) override
	{
		const int vertexAmount = stoi(parameters[1]);

		isOriented = parameters[2] == "1";
		isWeighted = parameters[3] == "1";

		adjacencyList = readList(stream, vertexAmount, isWeighted);
	}


	void addEdge(int from, int to, int weight) override
	{
		if (adjacencyList->at(from)->count(to) == 0 && (adjacencyList->at(to)->count(from) == 0 || isOriented))
		{
			++edgesAmount;
			++(*vertexDegrees)[from];
			(*(*adjacencyList)[from])[to] = new Edge(from, to, weight);
		}
	}


	void removeEdge(int from, int to)  override
	{
		--edgesAmount;
		removeEdgeSimplex(adjacencyList, from, to);
		if (!isOriented)
			removeEdgeSimplex(adjacencyList, to, from);
	}


	int changeEdge(int from, int to, int newWeight)  override
	{
		const int oldWeight = (*(*adjacencyList)[from])[to]->cap();
		(*(*adjacencyList)[from])[to]->setCap(newWeight);
		return oldWeight;
	}


	AdjMatrixGraphRepresentation* transformToAdjMatrix()  override;


	AdjListGraphRepresentation* transformToAdjList()  override { return this; }


	EdgeListGraphRepresentation* transformToListOfEdges()  override;


	void writeGraph(string fileName) override
	{
		ofstream file(fileName);
		file << "L " << adjacencyList->size() << endl;
		file << (isOriented ? "1 " : "0 ") << (isWeighted ? "1" : "0") << endl;

		writeList(file, adjacencyList, isWeighted);

		file.close();
	}


	int calcWeight() override {
		return 0;
	}
	int calcFlow() override {
		return 0;
	}


	AdjacencyList copyListWithDuplicates(bool ignoreOrientation) {
		AdjacencyList adjList = new vector<map<int, Edge*>*>(adjacencyList->size());
		for (size_t i = 0; i < adjacencyList->size(); ++i)
			(*adjList)[i] = new map<int, Edge*>();

		for (size_t i = 0; i < adjacencyList->size(); ++i)
			for (auto it = adjacencyList->at(i)->begin(); it != adjacencyList->at(i)->end(); ++it) {
				(*(*adjList)[i])[it->second->w()] = new Edge(i, it->second->w(), it->second->cap(), false);
				if (!isOriented || ignoreOrientation)
					(*(*adjList)[it->second->w()])[i] = new Edge(it->second->w(), i, it->second->cap(), true);
			}

		return adjList;
	}

	friend Graph;
};

class EdgeListGraphRepresentation : GraphRepresentation
{
	int vertexAmount;
	EdgeList edgeList = nullptr;
	EdgeList readList(istream& stream, int edgeAmount)
	{
		const EdgeList edgeList = new map<pair<int, int>, Edge*>();
		vertexDegrees = new map<int, int>();
		for (size_t i = 0; i < vertexAmount; ++i)
			(*vertexDegrees)[i] = 0;
		for (size_t i = 0; i < edgeAmount; ++i)
		{
			string line;
			getline(stream, line);
			vector<string> edgeParameters = StringParser::split(line, ' ');
			const int weight = edgeParameters.size() > 2 ? stoi(edgeParameters[2]) : 1;
			const int from = stoi(edgeParameters[0]) - 1;
			const int to = stoi(edgeParameters[1]) - 1;
			if (edgeList->count(make_pair(from, to)) == 0 && (edgeList->count(make_pair(to, from)) == 0 || isOriented))
			{
				(*edgeList)[make_pair(from, to)] = new Edge(from, to, weight);
				++(*vertexDegrees)[from];
				if (!isOriented)
					++(*vertexDegrees)[to];
			}
		}
		edgesAmount = edgeList->size();
		return edgeList;
	}
	static void writeList(ostream& stream, EdgeList edgeList, bool isWeighted)
	{
		for (auto it = edgeList->begin(); it != edgeList->end(); ++it)
		{
			stream << it->second->v() + 1 << " " << it->second->w() + 1;
			if (isWeighted)
				stream << " " << it->second->cap();
			if (next(it) != edgeList->end())
				stream << endl;
		}
	}
	void copy(const EdgeListGraphRepresentation& graph) {
		this->isOriented = graph.isOriented;
		this->isWeighted = graph.isWeighted;
		this->vertexAmount = graph.vertexAmount;
		this->edgesAmount = graph.edgesAmount;
		this->vertexDegrees = new map<int, int>(*graph.vertexDegrees);
		this->edgeList = new map<pair<int, int>, Edge*>(*(graph.edgeList));
	}
	int getEdgesAmount() {
		return edgeList->size();
	}
	DSU* buildDSU() {
		DSU* dsu = new DSU(vertexAmount);
		for (int i = 0; i < vertexAmount; ++i)
			dsu->makeSet(i);
		for (auto it = edgeList->begin(); it != edgeList->end(); ++it) {
			dsu->unite(it->second->v(), it->second->w());
		}
		return dsu;
	}

public:
	EdgeListGraphRepresentation() = default;
	EdgeListGraphRepresentation(bool isOriented, bool isWeighted, int vertexAmount, int edgesAmount)
	{
		this->vertexDegrees = new map<int, int>();
		this->edgeList = new map<pair<int, int>, Edge*>();
		this->isOriented = isOriented;
		this->isWeighted = isWeighted;
		this->vertexAmount = vertexAmount;
		this->edgesAmount = edgesAmount;
	}
	EdgeListGraphRepresentation(const EdgeListGraphRepresentation& graph) {
		copy(graph);
	}
	~EdgeListGraphRepresentation()
	{
		for (auto it = edgeList->begin(); it != edgeList->end(); ++it)
			delete(it->second);
		delete vertexDegrees;
		delete edgeList;
	}
	EdgeListGraphRepresentation& operator=(const EdgeListGraphRepresentation& graph) {
		copy(graph);
		return *this;
	}

	bool isAdjMatrixGraph() override { return false; }
	bool isAdjListGraph() override { return false; }
	bool isEdgeListGraph() override { return true; }

	void readGraph(istream& stream, vector<string> parameters) override
	{
		vertexAmount = stoi(parameters[1]);
		int edgeAmount = stoi(parameters[2]);

		isOriented = parameters[3] == "1";
		isWeighted = parameters[4] == "1";

		edgeList = readList(stream, edgeAmount);
	}


	void addEdge(int from, int to, int weight) override
	{
		if (edgeList->count(make_pair(from, to)) == 0 && (edgeList->count(make_pair(to, from)) == 0 || isOriented))
		{
			++edgesAmount;
			(*edgeList)[make_pair(from, to)] = new Edge(from, to, weight);
			++(*vertexDegrees)[from];
			++(*vertexDegrees)[to];
		}
	}


	void removeEdge(int from, int to)  override
	{
		--edgesAmount;
		auto it = edgeList->begin();
		bool edgeFound = false;
		for (; it != edgeList->end(); ++it) {
			if (it->second->v() == from && it->second->w() == to) {
				edgeFound = true;
				break;
			}
		}
		if (edgeFound) {
			delete(it->second);
			edgeList->erase(it);
			--(*vertexDegrees)[from];
			--(*vertexDegrees)[to];
		}
	}


	int changeEdge(int from, int to, int newWeight)  override
	{
		const pair<int, int> edge = make_pair(from, to);
		const int oldWeight = (*edgeList)[edge]->cap();
		(*edgeList)[edge]->setCap(newWeight);
		return oldWeight;
	}


	AdjMatrixGraphRepresentation* transformToAdjMatrix()  override;


	AdjListGraphRepresentation* transformToAdjList()  override;


	EdgeListGraphRepresentation* transformToListOfEdges()  override { return this; }


	void writeGraph(string fileName) override
	{
		ofstream file(fileName);
		file << "E " << vertexAmount << " " << edgeList->size() << endl;
		file << (isOriented ? "1 " : "0 ") << (isWeighted ? "1" : "0") << endl;

		writeList(file, edgeList, isWeighted);

		file.close();
	}

	friend Graph;

	int calcWeight() override {
		int result = 0;
		for (auto it = edgeList->begin(); it != edgeList->end(); ++it) {
			result += it->second->cap();
		}
		return result;
	}
	int calcFlow() override {
		int result = 0;
		for (auto it = edgeList->begin(); it != edgeList->end(); ++it) {
			if (it->second->flow() > 0)
				result += it->second->flow();
		}
		return result;
	}

	friend AdjListGraphRepresentation;
};

AdjListGraphRepresentation* AdjMatrixGraphRepresentation::transformToAdjList()
{
	AdjListGraphRepresentation* adjListGraph = new AdjListGraphRepresentation(adjacencyMatrix->size(), 0, isOriented, isWeighted);

	for (size_t i = 0; i < adjacencyMatrix->size(); ++i)
	{
		for (size_t j = 0; j < adjacencyMatrix->size(); ++j)
		{
			if ((*(*adjacencyMatrix)[i])[j] != nullptr)
			{
				int weight = (*(*adjacencyMatrix)[i])[j]->cap();
				if (weight > 0)
					adjListGraph->addEdge(i, j, weight);
			}
		}
	}

	return adjListGraph;
}

EdgeListGraphRepresentation* AdjMatrixGraphRepresentation::transformToListOfEdges()
{
	EdgeListGraphRepresentation* edgeListGraph = new EdgeListGraphRepresentation(isOriented, isWeighted, adjacencyMatrix->size(), 0);

	for (size_t i = 0; i < adjacencyMatrix->size(); ++i)
		for (size_t j = 0; j < adjacencyMatrix->size(); ++j)
		{
			if ((*(*adjacencyMatrix)[i])[j] != nullptr)
			{
				int weight = (*(*adjacencyMatrix)[i])[j]->cap();
				if (weight > 0)
					edgeListGraph->addEdge(i, j, weight);
			}
		}

	return edgeListGraph;
}

AdjMatrixGraphRepresentation* AdjListGraphRepresentation::transformToAdjMatrix()
{
	AdjMatrixGraphRepresentation* adjMatrixGraph = new AdjMatrixGraphRepresentation(adjacencyList->size(), 0, isOriented, isWeighted);

	for (size_t i = 0; i < adjacencyList->size(); ++i)
	{
		map<int, Edge*>* currentVertexMap = (*adjacencyList)[i];
		for (auto edge = currentVertexMap->begin(); edge != currentVertexMap->end(); ++edge) {
			adjMatrixGraph->addEdge(i, edge->first, edge->second->cap());
			adjMatrixGraph->adjacencyMatrix->at(i)->at(edge->first)->addflowRto(edge->first, edge->second->flow());
		}
	}

	return adjMatrixGraph;
}

EdgeListGraphRepresentation* AdjListGraphRepresentation::transformToListOfEdges()
{
	EdgeListGraphRepresentation* edgeListGraph = new EdgeListGraphRepresentation(isOriented, isWeighted, adjacencyList->size(), 0);

	for (size_t i = 0; i < adjacencyList->size(); ++i)
	{
		map<int, Edge*>* currentVertexMap = (*adjacencyList)[i];
		for (auto edge = currentVertexMap->begin(); edge != currentVertexMap->end(); ++edge) {
			edgeListGraph->addEdge(i, edge->first, edge->second->cap());
			edgeListGraph->edgeList->at(make_pair(i, edge->first))->addflowRto(edge->first, edge->second->flow());
		}
	}

	return edgeListGraph;
}

AdjMatrixGraphRepresentation* EdgeListGraphRepresentation::transformToAdjMatrix()
{
	AdjMatrixGraphRepresentation* adjMatrixGraph = new AdjMatrixGraphRepresentation(vertexAmount, 0, isOriented, isWeighted);

	for (auto it = edgeList->begin(); it != edgeList->end(); ++it)
		adjMatrixGraph->addEdge(it->first.first, it->first.second, it->second->cap());

	return adjMatrixGraph;
}

AdjListGraphRepresentation* EdgeListGraphRepresentation::transformToAdjList()
{
	AdjListGraphRepresentation* adjListGraph = new AdjListGraphRepresentation(vertexAmount, 0, isOriented, isWeighted);

	for (auto it = edgeList->begin(); it != edgeList->end(); ++it)
		adjListGraph->addEdge(it->first.first, it->first.second, it->second->cap());

	return adjListGraph;
}

class Graph
{
	GraphRepresentation* graphRepresentation = nullptr;
	static vector<string> getGraphParameters(istream &stream)
	{
		vector<string> parameters;
		string line;
		for (size_t i = 0; i < 2; i++)
		{
			getline(stream, line);
			for (const string token : StringParser::split(line, ' '))
				parameters.push_back(token);
		}
		return parameters;
	}
	static GraphRepresentation* decideRepresentationType(string code)
	{
		if (code == "C")
			return (GraphRepresentation*) new AdjMatrixGraphRepresentation();
		if (code == "L")
			return (GraphRepresentation*) new AdjListGraphRepresentation();
		if (code == "E")
			return (GraphRepresentation*) new EdgeListGraphRepresentation();
		return nullptr;
	}
	void copy(const Graph& graph) {
		GraphRepresentation* newGraph = nullptr;
		if (graph.graphRepresentation->isAdjListGraph())
			newGraph = (AdjListGraphRepresentation*)graph.graphRepresentation;
		else if (graph.graphRepresentation->isAdjMatrixGraph())
			newGraph = (AdjMatrixGraphRepresentation*)graph.graphRepresentation;
		else if (graph.graphRepresentation->isEdgeListGraph())
			newGraph = (EdgeListGraphRepresentation*)graph.graphRepresentation;
		delete graphRepresentation;
		graphRepresentation = newGraph;
	}
	MarkedEdges initMarkedEdges(AdjacencyList adjList) {
		MarkedEdges markedEdges = new vector<map<int, bool>*>(adjList->size());
		for (int i = 0; i < markedEdges->size(); ++i) {
			(*markedEdges)[i] = new map<int, bool>();
			for (auto it = adjList->at(i)->begin(); it != adjList->at(i)->end(); ++it)
				(*(*markedEdges)[i])[it->first] = false;
		}
		return markedEdges;
	}
	bool isBridge(AdjacencyList adjList, MarkedEdges traversed, int from, int to)
	{
		vector<int>* discoveryTime = new vector<int>(adjList->size(), -1);
		vector<int>* minTime = new vector<int>(adjList->size(), -1);
		vector<int>* parent = new vector<int>(adjList->size(), -1);

		int time = 0;
		return bridgeDFS(adjList, traversed, time, 0, from, to, discoveryTime, minTime, parent);
	}
	bool bridgeDFS(AdjacencyList adjList, MarkedEdges traversed, int& time, int source, int u, int v,
		vector<int>* discoveryTime, vector<int>* minTime, vector<int>* parent)
	{
		discoveryTime->at(source) = minTime->at(source) = time++;
		for (auto it = adjList->at(source)->begin(); it != adjList->at(source)->end(); ++it)
		{
			int next = it->first;
			if (!traversed->at(source)->at(next))
			{
				if (discoveryTime->at(next) == -1)
				{
					parent->at(next) = source;
					if (bridgeDFS(adjList, traversed, time, next, u, v, discoveryTime, minTime, parent))
						return true;
					minTime->at(source) = min(minTime->at(source), minTime->at(source));
					if (minTime->at(next) > discoveryTime->at(source) && (source == u && next == v || source == v && next == u))
						return true;
				}
				else if (next != parent->at(source))
					minTime->at(source) = min(minTime->at(source), discoveryTime->at(next));
			}
		}
		return false;
	}
	bool dfsKuhn(int v, vector<bool> & used, AdjacencyList adjList, vector<int> & part) {
		if (used[v])  return false;
		used[v] = true;
		for (auto it = adjList->at(v)->begin(); it != adjList->at(v)->end(); ++it) {
			int to = it->first;
			if (part[to] == -1 || dfsKuhn(part[to], used, adjList, part)) {
				part[to] = v;
				return true;
			}
		}
		return false;
	}
	int dfsFordFalkerson(int v, int Cmin, int sink, vector<bool>* visited, AdjacencyList adjList) {
		if (v == sink)
			return Cmin;
		visited->at(v) = true;
		for (auto it = adjList->at(v)->begin(); it != adjList->at(v)->end(); ++it) {
			Edge* vw = it->second;
			if (!visited->at(vw->w()) && vw->flowLeft() != 0) {
				int delta = dfsFordFalkerson(vw->w(), min(Cmin, vw->flowLeft()), sink, visited, adjList);
				if (delta > 0) {
					vw->addflowRto(vw->w(), delta);
					Edge* wv = adjList->at(vw->w())->at(vw->v());
					wv->addflowRto(vw->w(), delta);
					return delta;
				}
			}
		}
		return 0;
	}
	bool bfsDinic(int source, int sink, AdjacencyList adjList, vector<int>* distances) {
		for (int i = 0; i < distances->size(); ++i)
			distances->at(i) = INT32_MAX;
		queue<int> q = queue<int>();
		distances->at(source) = 0;
		q.push(source);
		while (!q.empty()) {
			int v = q.front();
			q.pop();
			for (auto it = adjList->at(v)->begin(); it != adjList->at(v)->end(); ++it) {
				if (it->second->flowLeft() != 0 && distances->at(it->first) == INT32_MAX) {
					distances->at(it->first) = distances->at(v) + 1;
					q.push(it->first);
				}
			}
		}
		return distances->at(sink) != INT32_MAX;
	}
	int dfsDinic(int v, int Cmin, int sink, vector<map<int, Edge*>::iterator>* p, vector<int>* distances, AdjacencyList adjList) {
		if (v == sink || Cmin == 0)
			return Cmin;
		for (; p->at(v) != adjList->at(v)->end(); ++p->at(v)) {
			Edge* vw = p->at(v)->second;
			if (distances->at(vw->w()) == distances->at(v) + 1) {
				int delta = dfsDinic(vw->w(), min(Cmin, vw->flowLeft()), sink, p, distances, adjList);
				if (delta != 0) {
					vw->addflowRto(vw->w(), delta);
					Edge* wv = adjList->at(vw->w())->at(v);
					wv->addflowRto(vw->w(), delta);
					return delta;
				}
			}
		}
		return 0;
	}
	void clearAdjList(AdjacencyList adjList) {
		for (int i = 0; i < adjList->size(); ++i) {
			for (auto it = adjList->at(i)->begin(); it != adjList->at(i)->end(); ++it)
				delete(it->second);
			delete(adjList->at(i));
		}
		delete(adjList);
	}
	Graph* createFlowGraph(AdjacencyList adjList)
	{
		AdjListGraphRepresentation* newGraphRepresentation = new AdjListGraphRepresentation(adjList->size(), graphRepresentation->edgesAmount,
			graphRepresentation->isOriented, graphRepresentation->isWeighted);
		for (int i = 0; i < adjList->size(); ++i)
			for (auto it = adjList->at(i)->begin(); it != adjList->at(i)->end(); ++it)
				if (!it->second->isDuplicate())
					newGraphRepresentation->addEdge(i, it->first, it->second->flow());
		Graph* graph = new Graph();
		graph->graphRepresentation = newGraphRepresentation;
		return graph;
	}

public:
	Graph() {}
	Graph(int n) {
		graphRepresentation = (GraphRepresentation*) new EdgeListGraphRepresentation(false, true, n, 0);
	}
	Graph(const Graph& graph) {
		copy(graph);
	}
	~Graph()
	{
		delete graphRepresentation;
	};
	Graph& operator=(const Graph& graph) {
		copy(graph);
		return *this;
	}

	void readGraph(string fileName)
	{
		ifstream file(fileName);
		const vector<string> parameters = getGraphParameters(file);
		graphRepresentation = decideRepresentationType(parameters[0]);
		graphRepresentation->readGraph(file, parameters);
		file.close();
	}
	void addEdge(int from, int to, int weight) const
	{
		graphRepresentation->addEdge(from, to, weight);
	}
	void removeEdge(int from, int to) const
	{
		graphRepresentation->removeEdge(from, to);
	}
	int changeEdge(int from, int to, int newWeight) const
	{
		return graphRepresentation->changeEdge(from, to, newWeight);
	}
	void transformToAdjList()
	{
		if (!graphRepresentation->isAdjListGraph()) {
			GraphRepresentation* adjListGraph = (GraphRepresentation*)graphRepresentation->transformToAdjList();
			delete graphRepresentation;
			graphRepresentation = adjListGraph;
		}
	}
	void transformToAdjMatrix()
	{
		if (!graphRepresentation->isAdjMatrixGraph()) {
			GraphRepresentation* adjMatrixGraph = (GraphRepresentation*)graphRepresentation->transformToAdjMatrix();
			delete graphRepresentation;
			graphRepresentation = adjMatrixGraph;
		}
	}
	void transformToListOfEdges()
	{
		if (!graphRepresentation->isEdgeListGraph()) {
			GraphRepresentation* edgeListGraph = (GraphRepresentation*)graphRepresentation->transformToListOfEdges();
			delete graphRepresentation;
			graphRepresentation = edgeListGraph;
		}
	}
	void writeGraph(string fileName) const
	{
		graphRepresentation->writeGraph(fileName);
	}
	//Алгоритм Прима
	Graph getSpaingTreePrima() {
		transformToAdjList();
		AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->adjacencyList;

		priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
		int src = 0;
		vector<int> distances(adjList->size(), INT32_MAX);
		vector<int> parent(adjList->size(), -1);
		vector<bool> inMST(adjList->size(), false);

		pq.push(make_pair(0, src));
		distances[src] = 0;

		while (!pq.empty()) {
			int u = pq.top().second;
			pq.pop();
			inMST[u] = true;
			for (auto it = adjList->at(u)->begin(); it != adjList->at(u)->end(); ++it) {
				int v = it->first;
				int weight = it->second->cap();
				if (inMST[v] == false && distances[v] > weight) {
					distances[v] = weight;
					pq.push(make_pair(distances[v], v));
					parent[v] = u;
				}
			}
		}

		EdgeListGraphRepresentation* newEdgeListGraph = new EdgeListGraphRepresentation(false, true, adjList->size(), 0);
		for (size_t i = 1; i < parent.size(); ++i)
			newEdgeListGraph->addEdge(parent[i], i, distances[i]);
		Graph* graph = new Graph();
		graph->graphRepresentation = newEdgeListGraph;
		return *graph;
	}
	//Крускал
	Graph getSpaingTreeKruscal() {
		transformToListOfEdges();
		EdgeListGraphRepresentation* edgeListGraph = (EdgeListGraphRepresentation*)graphRepresentation;

		auto cmp = [](pair<pair<int, int>, Edge*> a, pair<pair<int, int>, Edge*> b) {
			return a.second->cap() < b.second->cap();
		};
		auto edgeVector = vector<pair<pair<int, int>, Edge*>>(edgeListGraph->edgeList->begin(), edgeListGraph->edgeList->end());
		sort(edgeVector.begin(), edgeVector.end(), cmp);

		DSU dsu = DSU(edgeListGraph->vertexAmount);
		for (size_t i = 0; i < edgeListGraph->vertexAmount; ++i)
			dsu.makeSet(i);

		EdgeListGraphRepresentation* newEdgeListGraph = new EdgeListGraphRepresentation(false, true, edgeListGraph->vertexAmount, 0);
		for (size_t i = 0; i < edgeVector.size(); i++)
		{
			int from = edgeVector[i].first.first;
			int to = edgeVector[i].first.second;
			if (dsu.find(from) != dsu.find(to)) {
				newEdgeListGraph->addEdge(from, to, edgeVector[i].second->cap());
				dsu.unite(from, to);
			}
		}
		Graph* graph = new Graph();
		graph->graphRepresentation = newEdgeListGraph;
		return *graph;
	}
	//Борувка
	Graph getSpaingTreeBoruvka() {
		transformToListOfEdges();
		EdgeListGraphRepresentation* edgeListGraph = (EdgeListGraphRepresentation*)graphRepresentation;
		EdgeList edgeList = edgeListGraph->edgeList;
		int vertexAmount = edgeListGraph->vertexAmount;
		auto edgeVector = vector<pair<pair<int, int>, Edge*>>(edgeListGraph->edgeList->begin(), edgeListGraph->edgeList->end());

		DSU dsu = DSU(vertexAmount);
		for (size_t i = 0; i < edgeListGraph->vertexAmount; ++i)
			dsu.makeSet(i);

		EdgeListGraphRepresentation* newEdgeListGraph = new EdgeListGraphRepresentation(false, true, vertexAmount, 0);
		while (newEdgeListGraph->edgesAmount < vertexAmount - 1) {
			auto minEdges = map<int, int>();
			for (int i = 0; i < vertexAmount; ++i)
				minEdges[i] = -1;
			for (int i = 0; i < edgeVector.size(); ++i)
			{
				auto edge = edgeVector[i];
				int from = edge.first.first;
				int to = edge.first.second;
				int weight = edge.second->cap();
				int fromComponent = dsu.find(from);
				int toComponent = dsu.find(to);
				if (fromComponent != toComponent) {
					if (minEdges[fromComponent] == -1 || edgeVector[minEdges[fromComponent]].second->cap() > weight)
						minEdges[fromComponent] = i;
					if (minEdges[toComponent] == -1 || edgeVector[minEdges[toComponent]].second->cap() > weight)
						minEdges[toComponent] = i;
				}
			}
			for (int i = 0; i < minEdges.size(); i++) {
				if (minEdges[i] != -1) {
					pair<pair<int, int>, Edge*> edge = edgeVector[minEdges[i]];
					dsu.unite(edge.first.first, edge.first.second);
					newEdgeListGraph->addEdge(edge.first.first, edge.first.second, edge.second->cap());
				}
			}
		}

		Graph* graph = new Graph();
		graph->graphRepresentation = newEdgeListGraph;
		return *graph;
	}

	int checkEuler(bool &circleExist) {
		bool connectedComponentsCheck = graphRepresentation->connectedComponentsCheck();
		if (!connectedComponentsCheck) {
			circleExist = false;
			return 0;
		}
		return graphRepresentation->getBeginVertex(circleExist) + 1;
	}

	vector<int> getEuleranTourFleri() {
		bool circleExists;
		int currentVertex = checkEuler(circleExists) - 1;
		if (!circleExists && currentVertex == -1)
			return vector<int>();

		transformToAdjList();
		AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->copyListWithDuplicates(false);
		MarkedEdges isTraversed = initMarkedEdges(adjList);
		vector<int> path = vector<int>();

		path.push_back(currentVertex + 1);
		for (int i = 0; i < graphRepresentation->edgesAmount; ++i) {
			int nextVertex = 0;
			bool onlyBridges = true;
			for (auto it = adjList->at(currentVertex)->begin(); it != adjList->at(currentVertex)->end(); ++it) {
				if (!(*(*isTraversed)[currentVertex])[it->first]) {
					if (!isBridge(adjList, isTraversed, currentVertex, it->first)) {
						onlyBridges = false;
						nextVertex = it->first;
					}
					else if (onlyBridges)
					{
						nextVertex = it->first;
					}
				}
			}
			path.push_back(nextVertex + 1);
			(*(*isTraversed)[currentVertex])[nextVertex] = true;
			(*(*isTraversed)[nextVertex])[currentVertex] = true;
			currentVertex = nextVertex;
		}
		clearAdjList(adjList);
		return path;
	}

	vector<int> getEuleranTourEffective() {
		bool circleExists;
		int currentVertex = checkEuler(circleExists) - 1;
		if (!circleExists && currentVertex == -1)
			return vector<int>();

		transformToAdjList();
		AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->copyListWithDuplicates(false);
		MarkedEdges isTraversed = initMarkedEdges(adjList);
		vector<int> path = vector<int>();

		stack<int> s = stack<int>();
		s.push(currentVertex);
		while (!s.empty())
		{
			int w = s.top();
			for (auto it = adjList->at(w)->begin(); it != adjList->at(w)->end(); ++it) {
				if (!isTraversed->at(w)->at(it->first)) {
					s.push(it->first);
					(*isTraversed->at(w))[it->first] = true;
					(*isTraversed->at(it->first))[w] = true;
					break;
				}
			}
			if (w == s.top()) {
				s.pop();
				path.push_back(w + 1);
			}
		}
		clearAdjList(adjList);
		return path;
	}

	int checkBipart(vector<char> &marks) {
		transformToAdjList();
		AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->copyListWithDuplicates(false);

		marks = vector<char>(adjList->size(), -1);

		for (int i = 0; i < adjList->size(); ++i) {
			if (marks[i] == -1) {
				queue<int> q = queue<int>();
				q.push(i);
				marks[i] = 'A';
				while (!q.empty()) {
					int v = q.front();
					q.pop();
					for (auto it = adjList->at(v)->begin(); it != adjList->at(v)->end(); ++it) {
						int to = it->first;
						if (marks[to] == -1) {
							if (marks[v] == 'A')
								marks[to] = 'B';
							else
								marks[to] = 'A';
							q.push(to);
						}
						else if (marks[to] == marks[v]) {
							return 0;
						}
					}
				}
			}
		}
		clearAdjList(adjList);
		return 1;
	}

	vector<pair<int, int>> getMaximumMatchingBipart() {
		vector<char> marks = vector<char>();
		int isBipart = checkBipart(marks);
		if (isBipart) {
			vector<pair<int, int>> result = vector<pair<int, int>>();
			AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->copyListWithDuplicates(false);
			vector<int> part = vector<int>(adjList->size(), -1);

			for (size_t i = 0; i < adjList->size(); i++)
				if (marks[i] == 'A') {
					vector<bool> used = vector<bool>(adjList->size(), false);
					dfsKuhn(i, used, adjList, part);
				}
			for (size_t i = 0; i < adjList->size(); i++)
				if (marks[i] == 'B' && part[i] != -1)
					result.push_back(make_pair(part[i] + 1, i + 1));
			clearAdjList(adjList);
			return result;
		}
		return vector<pair<int, int>>();
	}

	Graph flowFordFulkerson(int source, int sink) {
		transformToAdjList();
		AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->copyListWithDuplicates(true);
		while (true) {
			vector<bool>* visited = new vector<bool>(adjList->size(), false);
			int flow = dfsFordFalkerson(source - 1, INT32_MAX, sink - 1, visited, adjList);
			delete(visited);
			if (flow == 0)
				break;
		}

		Graph* graph = createFlowGraph(adjList);
		clearAdjList(adjList);
		return *graph;
	}
	Graph flowDinitz(int source, int sink) {
		transformToAdjList();
		AdjacencyList adjList = ((AdjListGraphRepresentation*)graphRepresentation)->copyListWithDuplicates(true);
		vector<int>* distances = new vector<int>(adjList->size(), 0);
		while (bfsDinic(source - 1, sink - 1, adjList, distances)) {
			vector<map<int, Edge*>::iterator>* p = new vector<map<int, Edge*>::iterator>(adjList->size());
			for (int i = 0; i < adjList->size(); ++i)
				p->at(i) = adjList->at(i)->begin();
			while (true) {
				int flow = dfsDinic(source - 1, INT32_MAX, sink - 1, p, distances, adjList);
				if (flow == 0)
					break;
			}
			delete(p);
		}
		delete(distances);
		Graph* graph = createFlowGraph(adjList);
		clearAdjList(adjList);
		return *graph;	}

	int calcWeight() {
		return graphRepresentation->calcWeight();
	}
	int calcFlow() {
		return graphRepresentation->calcFlow();
	}
};