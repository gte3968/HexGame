//Daniel Bai 8/13/2020
//Initial version of Hex game, using Dijkstra shortest path algo to check who wins.

//Daniel Bai 07/31/2020
//Dijkstra's shortest path algorithm
//Implemented using connectivity matrix
//For undirected graph (edge matrix symmetric)
//For ease of manual checking, the weight type was set to int
//Graph is randomly generated using 4 parameters: number of nodes, edge density, min and max of edge weights
//The algorithm stops as soon as the destination is found, thus only have nodes between source and destination nodes
//It can be modified to continue running to generate the entire shortest path tree from a given source node

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <climits>
#include <queue>
#include <tuple>  //key methods: make_tuple, get<i>
#include <cstddef>
#include <chrono>
#include <random>
using namespace std;

ostream& operator<< (ostream& out, vector<int> vec){
  for(int i=0; i<vec.size(); i++)
    cout << vec[i] << ", ";
  cout << endl;
  return out;
}


//function returning a random number between 0 and 1
double prob()
{ //must do seed outside of this func and only once, else it's not random!
  return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

//return index of Min element of an array
int min_index(int* cost, int size)
{
  int min=INT_MAX;
  int minIdx;
  for(int i=0; i<size; i++){
    if(cost[i]<min)
      min = cost[i];
      minIdx = i;
  }
  return minIdx;
}

class graph{
public:
 graph(int V, double density, int dLow, int dHigh);
 ~graph();
 int get_num_node(){return V;}
 double get_density(){return density;}
 int get_dLow(){return dLow;}
 int get_dHigh(){return dHigh;}
 int get_num_edge();
 bool adjacent(int x, int y){if(edge[x][y]<get_dLow()) return true; else return false;}  //true if there is an edge between nodes x and y, false otherwise.
 int get_edge_value(int x, int y){return edge[x][y];}
 void set_edge_value(int x, int y, int val){edge[x][y] = edge[y][x] = val;}
 void add_edge(int x, int y, int weight);
 void delete_edge(int x, int y)
 {
   if(edge[x][y] != -1) edge[x][y] = edge[y][x] = -1;
   else cout << "No eddge exists between node " << x << " and " << y << endl;
 }  //remove an edge between nodes x and y
 vector<int> neighbors(int x);  //list of neighbors of node x
 void showAdjMatrix();
 int get_node_value(int x){return keys[x];}  //get key value of a node
 void set_node_value(int x, int a){keys[x] = a;}  //set key value of a node

protected:
 int V;  //num nodes, numbered 0 to V-1
 int E;  //num edges
 double density;  //edge density, range 0.2 and 0.4 for this assignment
 int dLow, dHigh;  //edge distance range 1 to 10 for this assignment
 int** edge;  //edge weight matrix, 2D array size VxV, thus 2 stars
 int* keys;  //key values for each node, array of size V, assuming nodes are labeled as integers starting 0
};

//graph constructor
graph::graph(int V=0, double density=0.0, int dLow=1, int dHigh=10):
   V(V), density(density), dLow(dLow),dHigh(dHigh)
{
 //edge density [0.0,1.0], 0.2 and 0.4 for this assignment. Edge distance range: low, high [1, 10] for this assignment)

 const int size = V;
 if(size == 0) return;  //no need to populate edge matrix

 keys = new int [size]; //init keys
 for(int i=0; i<size; i++){
   keys[i] = i;
 }

 edge = new int* [size];  //heap allocation
 for(int i=0; i<size; i++){
   edge[i] = new int [size];
 }
 //microseconds for rand seed. time() returns second, no good - if running many runs they all end up with same seconds thus same seed!
 uint64_t us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::
                  now().time_since_epoch()).count();
 srand(us);  //random seeed

 for(int i=0; i<size; i++){  //V: num Nodes of the graph, graph private member
   for(int j=0; j<size; j++){
      if(i==j)
       edge[i][j] = 0;
      else{
     //randomly initialize edges (0: no connection, else: random value between dLow and dHigh
        if(prob() < density){
          int tmp = static_cast<int>(dLow+prob()*(dHigh-dLow));
          edge[i][j] = edge[j][i] = tmp;
        }
        else
          edge[i][j] = edge [j][i] = 100000;  //-1 denotes no connection
       //cout << edge[i][j] << ",";
       //cout << get_edge_value(i,j) << ", ";
     }
   }
 }
}

int graph::get_num_edge(){
 const int size = V;
 int numEdge = 0;
 for(int i=0; i<size; i++){
   for(int j=0; j<size; j++){
     if(edge[i][j]>=get_dLow() && edge[i][j]<=get_dHigh()) ++numEdge;
   }
 }
 return numEdge/2;
}

void graph::add_edge(int x, int y, int weight){
 if(x==y){
   cout << "Can't add edge to self!" << endl;
   return;
 }
 if(edge[x][y] == -1)
   edge[x][y] = edge[y][x] = weight;
 else
   cout<<"Edge already exists between nodes " << x << " and " << y << "!" << endl;
 } //add an edge between nodes x and y

vector<int> graph::neighbors(int x){  //return all neighbors of node x in a vector
 vector<int> neighborList;
 const int size = V;
 for(int i=0; i<size; i++){
   if(x != i){  //excluding self
     if(edge[x][i]>=get_dLow() && edge[x][i]<=get_dHigh())
       neighborList.push_back(i);
   }
 }
 return neighborList;
}

void graph::showAdjMatrix(){
  const int size = V;
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      cout << get_edge_value(i,j) << ", ";
    }
    cout << endl;
  }
}


graph::~graph(){  //destructor
 const int size = V;
// for(int i=0; i<size; i++){
//   delete [] edge[i];  //delete only takes pointer
// }
//delete [] edge;
}


//Below is class PriorityQueue
typedef pair<int, int> cost_node;  //cost_node is a type name of a pair, 1st is cost, 2nd is node ID
//typedef priority_queue<cost_node, vector<cost_node>, greater<cost_node>> cost_node_PQ;

class PriorityQueue{  //augmented from STL priority_queue, the built-in methods top(), push(), pop() still usable
  public:
    PriorityQueue(){}  //constructor
    int size();
    priority_queue<cost_node, vector<cost_node>, greater<cost_node>>& getPQ(){return pq;}  //return reference so inplace mutation can be done to pq
    bool contains(int node);
    void updatePriority(cost_node pair);  //update the priority queue if a node's cost is changed
    void push(cost_node pair){pq.push(pair);}
    void pop(){pq.pop();}
    bool empty(){return pq.empty();}
    cost_node top(){return pq.top();}
    void show();

  private:
    priority_queue<cost_node, vector<cost_node>, greater<cost_node>> pq;
};

int PriorityQueue::size(){  //return size of pq while keeping pq intact
  int s=0;
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;

  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    tmpPairArr.push_back(make_pair(fst, snd));
    pq.pop();
  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);

  return tmpPairArr.size();
}

//return true if the pair with the specified node is in the pq
bool PriorityQueue::contains(int node){
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;
  bool flag = false;

  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    tmpPairArr.push_back(make_pair(fst, snd));
    pq.pop();
    if(snd == node){
     flag=true;
     break;
    }
  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);

  return flag;
}

void PriorityQueue::updatePriority(cost_node pair){  //update the priority queue with the fed pair (cost, node), i.e. update cost for the node
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;

  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    if(snd == pair.second){
      fst = pair.first;
      snd = pair.second;
     }
    tmpPairArr.push_back(make_pair(fst, snd));
    pq.pop();
    if(snd == pair.second)
      break;

  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);
}

void PriorityQueue::show(){
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;

  if(pq.empty()){
    cout << "The priority queue is empty!" << endl;
    return;
  }
  cout << "PriorityQueue:" << endl;
  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    cout << "(" << fst << ", " << snd << ")" << endl;
    pq.pop();
    tmpPairArr.push_back(make_pair(fst, snd));
  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);
}

/************************************************************
class for Dijkstra shortest path
**************************************************************/
class ShortestPath{  //find shortest path for a given undirected graph g
  public:
    ShortestPath(graph& g, int u, int w, bool pathFound):g(g), u(u), w(w), pathFound(pathFound){}
    graph& getGraph(){return g;}
    int getSrcNode(){return u;}
    int getDestNode(){return w;}
    vector<int> path(int u, int w);
    vector<int> vertices();  //call graph g by reference (when calling, just use g, not &g)
    void showPath(const vector<int>& path);
    void setTotCost(int cost){totCost = cost;}
    int getTotCost(){return totCost;}
    void showCost();
    int getNumSteps(){return numSteps;}
    void setNumSteps(int nSteps){numSteps = nSteps;}
    bool isPathFound(){return pathFound;}
    void setPathFound(){pathFound = true;}

  private:
    graph g;
    int u, w;
    int totCost; //cost from u to w shortest path
    int numSteps; //total steps from u to w shortest path
    bool pathFound = false;

};

//returns the list of vertices of a graph g
vector<int> ShortestPath::vertices(){
  vector<int> nodes;
  for(int i=0; i<getGraph().get_num_node(); i++){
    nodes.push_back(i);
  }
  return nodes;
}

//find shortest path between source node u and destination node w and returns the sequence of vertices representing shorest path u-v1-v2-…-vn-w.
vector<int> ShortestPath::path(int srcNode, int destNode){  //u: source
  const int numV = getGraph().get_num_node();  //num vertices
  PriorityQueue openSet ;  //for nodes being considered but not finalized yet. See PriorityQueue class
  vector<cost_node> closeSet;
  int closeSetArr[numV];  //flag for nodes in closeSet
  int costList[numV];  //cost of each node to source
  int prev[numV];  //prev of each node along shortest path
  for(int i=0; i<numV; i++){  //init
    costList[i] = 100000;
    closeSetArr[i] = false;
    prev[i] = -1;  //initial no parents. When building the path later, need to ensure idx of prev != -1
    if(i == srcNode)
      costList[i]=0;
  }

//populate openSet with source node
  int currNode = srcNode;
  int tmpCost = costList[currNode];
  cost_node currPair = make_pair(tmpCost, currNode);
  openSet.push(currPair);

  while(!openSet.empty()){
    //openSet.show();
    cost_node topPair = openSet.top();
    tmpCost = topPair.first;
    currNode = topPair.second;
    closeSet.push_back(topPair);

    closeSetArr[currNode] = true;  //flag true if a node is put in closeSet
    openSet.pop();
    if(currNode == destNode) break;  //reached w, done

    vector<int> nnList = getGraph().neighbors(currNode);  //get all nearest neighbors of the node that was just placed in closeSet
    for(int i=0; i<nnList.size(); i++){  //find shortest path from currNode to its NN
      int newNode = nnList[i];
      int parent = currNode;
      if(closeSetArr[newNode]) continue;  //if the node is already in closeSet skip it
      int currCost = getGraph().get_edge_value(currNode, newNode);
      cost_node newPair;

       if(tmpCost + currCost < costList[newNode]){
        costList[newNode] = tmpCost + currCost;  //relaxation
        prev[newNode] = parent;
        newPair = make_pair(costList[newNode], newNode);

        if(!(openSet.contains(newNode)))
          openSet.push(newPair);
        else
          openSet.updatePriority(newPair);
       }
    }
  }

  vector<int> returnedPath;  //for function return
  if(closeSet[closeSet.size()-1].second != destNode){  //if w was found it should be the last of closeSet
    //cout << "No path found between node " << srcNode << " and node " << destNode << "!" << endl;
    vector<int> emptyVec;
    return emptyVec;  //no path found
  }
  else{
    vector<int> path;
    currNode = destNode;
    path.push_back(currNode);
    while(prev[currNode]!=-1 && prev[currNode] <= numV){  //-1 indicates no parents
      path.push_back(prev[currNode]);
      currNode = prev[currNode];
    }
    //reverse order of path to src to dest
    for(int i=path.size()-1; i>=0; i--)
      returnedPath.push_back(path[i]);
    //save total cost
    setTotCost(costList[destNode]);
    setNumSteps(path.size()-1);
  }
  setPathFound();
  //cout << "pathFound = " << pathFound << endl;
  return returnedPath;
}

void ShortestPath::showPath(const vector<int>& path){
  cout << "The shortest path from node " << getSrcNode() << " to node " << getDestNode() << " is: " << endl;
  for(int i=0; i<path.size(); i++){
    cout << path[i];
    if(i!=path.size()-1)
      cout << "->";
  }
}

void ShortestPath::showCost(){
  int cost;
  cost = getTotCost();
  cout << "\nTotal cost from node " << getSrcNode() << " to node " << getDestNode() << ": " << cost << endl;

}

/*************************************************************************
Hex game
**************************************************************************/

#include <iostream>
#include <utility>  //pair
#include <cassert>
#include <algorithm>
using namespace std;
//enum class COLOR{BLUE, RED};

class Board{
  public:
    Board(int size);
    void drawHexBoard();
    pair<int, int> hexCoorToBoardCoor(int i, int j){return make_pair(2*i, 4*j+2*i);}
    void placeHex(int& i, int& j, char ch); //place a hex (represented by 'X' and 'O') at position (i,j), from (0,0) to (n-1, n-1)
        //Note there is a translation from position of node (i,j) to its position in the grid 2D array done by hexCoorToBoardCoor()
    bool isTaken(int i, int j);  //whether a node at position (i,j) has been taken

  private:
    int n; //nxn hex grids
    char** grid;  //grid is a 2D Vector of size 6*(n-1)+1 (per row) by 2*(n-1)+1 (# rows)
};

//Initialize the grid 2D array for the board
Board::Board(int size){ //initialize the grid 2D array
   n = size;
   const int row = 2*(n-1)+1;  //num rows
   const int col = 6*(n-1)+1;  //num chars per row
   grid = new char* [row];  //2D array allocation
   for(int i=0; i<row; ++i)
     grid[i] = new char [col];

   for(int i=0; i<row; ++i)  //init all chars to space
     for(int j=0; j<col; ++j)
       grid[i][j] = ' ';


   for(int i=0; i<row; ++i){
     for(int j=0; j<col; ++j){
       if(!(i%2)){  //rows with hex (dots, or nodes)
         for(int j=i; j<i+4*(n-1)+2; ++j){
           if((j-i)%4==0)
             grid[i][j] = '.';
           else if((j-i)%4==2)
             grid[i][j] = '-';
         }
       }
       else{  //rows with / and \ for edges
         for(int j=i; j<i+4*(n-1)+2; ++j){
           if((j-i)%4==0)
             grid[i][j] = '\\';
           else if((j-i)%4==2)
             grid[i][j] = '/';
         }
       }
     }
  }
}

//draw the board
void Board::drawHexBoard(){
  const int row = 2*(n-1)+1;  //num rows
  const int col = 6*(n-1)+1;  //num chars per row
  cout << endl;
  for(int i=0; i<row; ++i){
    for(int j=0; j<col; ++j){
      cout << grid[i][j];
    }
    cout << endl;
  }
}

//place a hex (represented by 'X' and 'O') at position (i,j), from (0,0) to (n-1, n-1)
void Board::placeHex(int& i, int& j, char ch){
  pair<int, int> gridCoor;
  int iBoard, jBoard;
  gridCoor = hexCoorToBoardCoor(i,j);  //convert hex coor to position on Board
  iBoard = gridCoor.first, jBoard=gridCoor.second;
  int ii=i, jj=j;
  while(isTaken(ii, jj)){
    cout << "Illegal move. Try again..." << endl;
    cin >> ii >> jj;
  }
  gridCoor = hexCoorToBoardCoor(ii,jj);  //convert hex coor to position on Board
  iBoard = gridCoor.first, jBoard=gridCoor.second;
  grid[iBoard][jBoard] = ch;
}

//whether a node at position (i,j) has been taken
bool Board::isTaken(int i, int j)
{
  pair<int, int> gridCoor;
  gridCoor = hexCoorToBoardCoor(i,j);
  int iBoard = gridCoor.first, jBoard=gridCoor.second;
  if(grid[iBoard][jBoard] != '.')
    return true;
  else
    return false;
}


void draw_Hex_Board(int n) //n is size of borad, for 7x7 hex board, size = 7;
{
  for(int i=0; i<2*n-1; i++){
    for(int k=0; k<i; ++k)
      cout << " ";
    if(!(i%2)){
      for(int j=0; j<n; ++j){
        if(j==n-1)
          cout << ".";
        else
          cout << ". - ";
      }
      cout << "\n";
    }
    else{
     for(int j=0; j<n; ++j){
       if(j==n-1)
         cout << "\\";
       else
         cout << "\\ / ";
     }
     cout << "\n";
    }
  }
}

/********************************************************
inherit graph class from Part A Dijkstra homework for here
*********************************************************/
class HexGraph: public graph{
  public:
    HexGraph(int s): graph(s*s, 0, 1, 1), size(s){};  //graph not connected at all to begin width, each time a node is filled, update edges
    int HexCoorToNode(int i, int j){return i*size+j;}  //convert hex coordinate to node #
    pair<int, int> NodeToHexCoor(int node){int i=node/size; int j=node%size; return make_pair(i,j);}
    vector<int> NeighborHexList(int i, int j);  //i,j: hex coordinate
    bool isGameOver();
    int GetSize(){return size;}
  private:
    int size;  //if size =7, it's a 7x7 hexBoard
};

//hex neighbor nodes
vector<int> HexGraph::NeighborHexList(int i, int j)
{
  int nnNode;
  vector<int> nnList;
  int nextI, nextJ;
  int flag; // 1: (0,0), 2:(0, size-1), 3:(size-1, 0), 4:(size-1, size-1), 5:(0,y), 6:(size-1,y), 7:(x,0), 8:(x, size-1), 9:(x,y)
  if(i==0){  //1-4; corner
    if(j==0) flag=1;
    else if(j==size-1) flag=2;
    else flag=5;  //5,6: x edge
  }
  else if(i==size-1){
    if(j==0) flag=3;
    else if(j==size-1) flag=4;
    else flag=6;
  }
  else if(j==0){  //7,8: y edge
    flag=7;
  }
  else if(j==size-1){
    flag=8;
  }
  else{  //internal
    flag=9;
  }

  switch(flag)
  {
    case 1:{
      nnNode = HexCoorToNode(i+1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j+1); nnList.push_back(nnNode);
      return nnList;
    }
    case 2:{
      nnNode = HexCoorToNode(i+1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j-1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i+1,j-1); nnList.push_back(nnNode);
      return nnList;
    }
    case 3:{
      nnNode = HexCoorToNode(i-1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j+1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j+1); nnList.push_back(nnNode);
      return nnList;
    }
    case 4:{
      nnNode = HexCoorToNode(i-1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j-1); nnList.push_back(nnNode);
      return nnList;
    }
    case 5:{  //i=0
      nnNode = HexCoorToNode(i,j+1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j-1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i+1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i+1,j-1); nnList.push_back(nnNode);
      return nnList;
    }
    case 6:{  //i=size-1
      nnNode = HexCoorToNode(i,j+1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j-1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j+1); nnList.push_back(nnNode);
      return nnList;
    }
    case 7:{  //j==0
      nnNode = HexCoorToNode(i+1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j+1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j+1); nnList.push_back(nnNode);
      return nnList;
    }
    case 8:{  //j=size-1
      nnNode = HexCoorToNode(i+1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j-1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i+1,j-1); nnList.push_back(nnNode);
      return nnList;
    }
    case 9:{  //internal hex
      nnNode = HexCoorToNode(i-1,j); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i-1,j+1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j-1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i,j+1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i+1,j-1); nnList.push_back(nnNode);
      nnNode = HexCoorToNode(i+1,j); nnList.push_back(nnNode);
      return nnList;
    }
  }
}

bool HexGraph::isGameOver(){
   //TODO detect last move player, then only check for that Player
   //For now check both players graphs using Dijkstra algo to see if it's connected edge to edges

   //check player 1 (top-down)


}

/********************************************************
class Player
*********************************************************/
class Player{
  public:
    friend class Board; //to access size of Board
    Player(int i): id(i){assert (id==1 || id==2 ||id==3); setSymbol();} //3 is the tmpP2 for MC
    void setSymbol(){if(id==1) symbol = 'X'; else if(id==2) symbol = 'O';}
    char getSymbol(){return symbol;}
    int getPlayerID(){return id;}
    pair<int, int> Move(HexGraph& g, Board& hBoard, Player opponent);  //make a move, returns <i,j> the position of the hex just placed
    void UpdateHexPlayed(int i, int j){hexPlayed.push_back(make_pair(i,j));}
    void UpdateNodePlayed(int i, int j, int size){nodePlayed.push_back(i*size+j);}
    void SetNodePlayed(vector<int> vec){nodePlayed = vec;}
    void UpdateEdgeList(HexGraph& g, int& i, int& j);
    vector<pair<int, int>> GetHexPlayed(){return hexPlayed;}
    vector<int> GetNodePlayed(){return nodePlayed; }
    void UpdateEdgeNodeList(HexGraph& g, int& i, int& j);
    vector<int> GetFirstEgeNodeList(){return firstEdgeNodeList;}
    vector<int> GetSecondEgeNodeList(){return secondEdgeNodeList;}
    bool IfWon(HexGraph& g);
    int runMC(HexGraph& h, vector<int> p1NodePlayed, vector<int> p2NodePlayed);

  private:
    int id;  //player 1 or 2
    char symbol;
    vector<pair<int, int>> hexPlayed;
    vector<int> nodePlayed;
    vector<int> firstEdgeNodeList, secondEdgeNodeList;  //track hex placed on the edge for checking if a player has won
};

//Make a move, including these steps:
//get the (i,j) and take input again if illegal move
//place the hex on the Board
//UpdateHexPlayed
//UpdateNodePlayed
//UpdateEdgeList
//returns the pair (i,j) just placed

pair<int, int> Player::Move(HexGraph& g, Board& hBoard, Player opponent)
{
  int i, j;
  if(id==1){
    cout << "Enter the position (i,j) to place your hex... " << endl;  //Player move
    int iP1, jP1;
    cin >> iP1 >> jP1;
    hBoard.placeHex(iP1, jP1, getSymbol());
    UpdateHexPlayed(iP1, jP1);
    UpdateNodePlayed(iP1,jP1, g.GetSize());
    UpdateEdgeList(g, iP1, jP1);
    i=iP1, j=jP1;
  }
  else if(id==2){
    int nextMoveNode = runMC(g, opponent.GetNodePlayed(), GetNodePlayed());
    pair<int, int> nextIJ = g.NodeToHexCoor(nextMoveNode);;
    int iP2, jP2;
    iP2 = nextIJ.first, jP2 = nextIJ.second;
//    while(true){  //computer's turn. if (i,j) is taken, try next pair. This version computer plays randomly
//      iP2 = rand()%g.GetSize();
//      jP2 = rand()%g.GetSize();
//      if(!hBoard.isTaken(iP2,jP2)){
        hBoard.placeHex(iP2, jP2, getSymbol());
        UpdateHexPlayed(iP2, jP2);
        UpdateNodePlayed(iP2, jP2, g.GetSize());
        UpdateEdgeList(g, iP2, jP2);
        i=iP2, j=jP2;
//      }
  }
  return(make_pair(i,j));
}

void Player::UpdateEdgeList(HexGraph& g, int& i, int& j)
{
  //for(auto it: this->GetHexPlayed()){cout<<"("<<it.first<< ","<<it.second<<") "<<endl;}
  int currNode = g.HexCoorToNode(i,j);
  vector<int> nnVec;
  nnVec = g.NeighborHexList(i,j);  //all neighbors of the hex just placed
  for(auto it:nnVec){  //update edges between the newly placed hex and already played hexes, if they are NN
      vector<int> tmpVec = GetNodePlayed();  //nodes placed by this player
      int tmpSize = tmpVec.size();
      for(int i=0; i<tmpSize; ++i){
        if(tmpVec[i] == it){ //one of the NN of newly placed hex has been placed previously
        g.set_edge_value(currNode, it, 1);  //undirected, so both ways
        g.set_edge_value(it, currNode, 1);
      }
    }
  }

}

//push (i,j) to the EdgeNodeList, the nodes on the edges of the board
void Player::UpdateEdgeNodeList(HexGraph& g, int& i, int& j)
{
  int flagIJ;
  if(id==1)  flagIJ=j;  //player1 go left-right, so check j for connection
  else if(id==2) flagIJ=i;  //player 2 go top-down, so check i for connection
  if(flagIJ==0){
    firstEdgeNodeList.push_back(g.HexCoorToNode(i,j));
  }
  if(flagIJ==g.GetSize()-1){
    secondEdgeNodeList.push_back(g.HexCoorToNode(i,j));
  }
}

bool Player::IfWon(HexGraph& g)
{
  for(auto ptr1=firstEdgeNodeList.begin(); ptr1<firstEdgeNodeList.end(); ++ptr1){  //any path between any node on the first edge and second edge
    for(auto ptr2=secondEdgeNodeList.begin(); ptr2<secondEdgeNodeList.end(); ++ptr2){
      int edge1Node = *ptr1;
      int edge2Node = *ptr2;
      ShortestPath sp1(g,edge1Node, edge2Node, false);
      sp1.path(edge1Node, edge2Node);
      if(sp1.isPathFound()){
        if(id==1){   //player
          cout << "You won!" << endl;
          return true;
        }
        else if(id==2){   //Computer
          cout << "You lost!" << endl;
          return true;
        }
      }
    }
  }
  return false;
}

/**********************************************************************************************************
//Run Monte Carlo for the computer given the nodes played so far by player and computer,
//Computer tries all available positions as next move candidates, then fill up the board randomly, and determines who won,
//run many times for each tentative position, determine the best move and return the node for the move
**********************************************************************************************************/
int Player::runMC(HexGraph& h, vector<int> p1NodePlayed, vector<int> p2NodePlayed)
{
  //remaining open nodes
  int size = h.GetSize();
  vector<int> emptyNodes(size*size, 1), availNovesVec;;  //1: empty, 0: taken
  for(auto it: p1NodePlayed) emptyNodes[it] = 0;  //initialize
  for(auto it: p2NodePlayed) emptyNodes[it] = 0;
  for(int i=0; i<size*size; ++i){
    if(emptyNodes[i]==1) availNovesVec.push_back(i);
  }
  cout << "available Nodes:";
  for(auto it: availNovesVec) cout << it << ",";
  cout << endl;

  srand(time(0));
  int numMCRun = 399;
  int bestMoveNode;

  vector<int> firstEdgeNodeList; secondEdgeNodeList;
  for(int i=0; i<size; ++i){ //define edge nodes for computer: going left to right
    firstEdgeNodeList.push_back(i);
    secondEdgeNodeList.push_back(size*(size-1)+i);
  }

  int numAvailHex, numAvailMoveP2; //number of available positions, and moves for computer after computer takes the tentative Move
  numAvailHex = size*size - p1NodePlayed.size() - p2NodePlayed.size() - 1;
  numAvailMoveP2 = numAvailHex / 2;
  cout << "numAvailMoveP2: " << numAvailMoveP2 << endl;

  priority_queue <pair<double, int>> winPercQueue;  //winPercQueue (winPerc, Node) for all the trials
  for(int trialNode=0; trialNode< emptyNodes.size(); ++trialNode){ //try all remaining available positions as next move candidates, run MC on each to get win%
    //HexGraph g(size);
    Player tmpP2(3);

    if(emptyNodes[trialNode]==0) continue;  //skip the node if already taken
    cout << "Trial Node: " << trialNode << ",";  //i is the next tentative move Node
    int winCount = 0;

    for(int j=0; j<numMCRun; ++j){  //MC for the current candidate move
      HexGraph g(size);
      tmpP2.SetNodePlayed(p2NodePlayed);  //copy current node played list of computer to the tmpP2 player
      //cout << "\nNodes played by computer so far: " << endl;
      //cout << tmpP2.GetNodePlayed() << endl;

  //TODO: somehow p2.Nodeplayed didn't get copied into tmpP2.nodesplayed
      for(auto it: tmpP2.GetNodePlayed()){  //update edges for the nodes played so far by computer
        //cout << it << ", ";
        pair<int, int> tmpNode = g.NodeToHexCoor(it);
        int tmpI = tmpNode.first, tmpJ = tmpNode.second;
        tmpP2.UpdateEdgeList(g, tmpI, tmpJ);
      }
      pair<int, int> tmpNode = g.NodeToHexCoor(trialNode);  //add trial node to tmpP2's played nodes list and update edges
      int tmpI = tmpNode.first, tmpJ = tmpNode.second;
      tmpP2.UpdateNodePlayed(tmpI, tmpJ, g.GetSize());
      tmpP2.UpdateHexPlayed(tmpI, tmpJ);
      tmpP2.UpdateEdgeList(g, tmpI, tmpJ);
      //g.showAdjMatrix();

      random_shuffle(availNovesVec.begin(), availNovesVec.end());  //shuffle remaining openings
      //for(auto it: availNovesVec) cout << it << " ";
      int counter =0;
      auto ptr = availNovesVec.begin();
      while(counter < numAvailMoveP2){  //randomly assign half of the opending positions to computer
        int theNode = *(ptr++);
        pair<int, int> tmpNode = g.NodeToHexCoor(theNode);
        if(theNode!=trialNode){  //if the Node is the trialNode, skip it, don't count up
          int tmpI = tmpNode.first, tmpJ = tmpNode.second;
          tmpP2.UpdateNodePlayed(tmpI, tmpJ, g.GetSize());
          tmpP2.UpdateHexPlayed(tmpI, tmpJ);
          tmpP2.UpdateEdgeList(g, tmpI, tmpJ);
          counter++;
        }
      }

      //cout << "tmpP2 nodes played: " << endl;
      //for(auto it: tmpP2.GetNodePlayed()) cout << it << ", ";
      //cout << endl;

      for(auto ptr1=firstEdgeNodeList.begin(); ptr1<firstEdgeNodeList.end(); ++ptr1){  //check if computer won for this MC run
        int flag=0;
        for(auto ptr2=secondEdgeNodeList.begin(); ptr2<secondEdgeNodeList.end(); ++ptr2){
          int edge1Node = *ptr1;
          int edge2Node = *ptr2;
          ShortestPath sp(g,edge1Node, edge2Node, false);
          //cout << "isPathFound before running path():" << sp.isPathFound() << endl;
          //cout << "edge1Node, edge2Node: " << edge1Node <<"," << edge2Node << endl;
          sp.path(edge1Node, edge2Node);
          if(sp.isPathFound()){
            winCount++;
            flag = 1;
            break;
           }
          }
        if(flag == 1) break;
       }
     }
    double winPerc = 1.0 * winCount / numMCRun;
    winPercQueue.push(make_pair(winPerc, trialNode));
    cout << "(win%, node): " << winPerc<< ", " << trialNode << endl;
    //cout << "winPerc, trialNode: " << winPerc << "," << trialNode << endl;
    //cout << "winPercQueue: " << winPercQueue.top().first << ", " << winPercQueue.top().second << endl;
    //cout << "winPercQueue.size(): " << winPercQueue.size() << endl;
   }

   //return the top winning % move
    if(!winPercQueue.empty()){
      pair<double, int> topMove = winPercQueue.top();
      cout << "topMove (win%, node): " << topMove.first << ", " << topMove.second << endl;
      bestMoveNode = topMove.second;
    }

 return bestMoveNode;
}

/****************************************************************
Driver code
*****************************************************************/

int main()
{
  int size = 7;
  cout << "Welcome to the Hex game." << endl << "You go first, connect left to right. Computer goes second, connect top down." <<endl;
  cout << "Enter the size n [3,11] for the n x n board:  " << endl;
  cin >> size;
  assert(size>=3 && size<=11);

  Board hexBoard(size);
  Player p1(1), p2(2);
  HexGraph g1(size);
  int i, j;
  srand(time(0));

  while(true){
    pair<int, int> ij = p1.Move(g1, hexBoard, p2);
    int iP1 = ij.first, jP1 = ij.second;
    //runMC(HexGraph& h, vector<int> p1NodePlayed, vector<int> p2NodePlayed);
    //hexBoard.drawHexBoard();
    p1.UpdateEdgeNodeList(g1, iP1, jP1);
    if(p1.IfWon(g1)) break;

    ij = p2.Move(g1, hexBoard, p1);
    hexBoard.drawHexBoard();
    int iP2 = ij.first, jP2 = ij.second;
    p2.UpdateEdgeNodeList(g1, iP2, jP2);
    if(p2.IfWon(g1)) break;
  }
  return 0;
}
