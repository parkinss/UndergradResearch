#include "SeedMatrix.hpp"

SeedMatrix::SeedMatrix(MeasureCombination * MC, double delta, bool isMaxHeap, std::unordered_set<ushort> & lex, std::unordered_set<ushort> & rex) :
  delta(delta), isMaxHeap(isMaxHeap),
  left_exclude(lex), right_exclude(rex)
{
  sims = MC->getAggregatedLocalSims();
  init_ptr = 0;
  init_column_vector();
}
SeedMatrix::~SeedMatrix(){
};

void SeedMatrix::init_column_vector(){
  std::cout << "init column vector" << std::endl;
  column_vector.reserve(sims.size());
  for(ushort i = 0; i < sims.size(); ++i){
    float max = -1;
    for(ushort j = 0; j < sims[i].size(); ++j){
      if(sims[i][j] > max){
	max = sims[i][j];
      }
    }
    seed_t curr;
    curr.left_node = i;
    curr.pairs = NULL;
    curr.top = max;
    curr.start = 0;
    curr.end = (uint)sims[i].size();
    /*
    seed_t curr = {
      i, 
      NULL,
      max,
      0, 
      (uint)sims[i].size()
    };
    */
    column_vector.push_back(curr);
  }
  //use reverse vectors to get a descending sort
  std::sort(column_vector.rbegin(), column_vector.rend());
  //this->debug();
}

std::pair<ushort,ushort> SeedMatrix::pop_uniform(){
  if(column_vector.size() <= 0){
    return std::make_pair(0,0);
    //throw QueueEmptyException();
  }

  std::pair<ushort,ushort> to_return;
  uint left_index;
  do {

  float tail = (column_vector[0].top - delta);

  /*
  seed_t curr = 
  
  std::cout << (curr.pairs ? "exist" : "NULL") << std::endl;
  std::cout << "(" << curr.getStart() << ", " 
	    << curr.getEnd() << ")" << std::endl;
  //don't run a search out here: the top row doesn't exist yet!
  //uint tmp = bin_search(*(curr.pairs), tail, curr.getStart(), curr.getEnd()); 
  //curr.end 
  */    
  uint i = init_ptr;
  for(; i < column_vector.size() && column_vector[i].top > tail; ++i){
    if(!column_vector[i].pairs){
      column_vector[i].pairs = create_node(column_vector[i].left_node);
      column_vector[i].end = bin_search(*(column_vector[i].pairs),
					tail,
					column_vector[i].getStart(),
					column_vector[i].getEnd());
    }
  }
  init_ptr = i;
  std::cout << "initialization done" << std::endl;
  uint choice_count = 0;
  for(uint i = 0; i < init_ptr; ++i){
    choice_count += (column_vector[i].end - column_vector[i].start);
  }
  std::uniform_int_distribution<int> uni(0, choice_count);
  std::random_device rd;
  std::mt19937 gen(rd());
  uint choice = uni(gen);
  std::cout << "pick random " << choice_count << " and " << choice << std::endl;
  uint counter = 0;
  while (choice > (column_vector[counter].end - column_vector[counter].start)) {
    choice -= (column_vector[counter].end - column_vector[counter].start);
    counter++;
  }
  //zero indexing (0, n-1)
  //std::uniform_int_distribution<int> uni(0, choice_count - 1);

  to_return = std::make_pair(column_vector[counter].left_node,
    column_vector[counter].pairs->at(column_vector[counter].start + choice).node);
  left_index = counter;

  column_vector[counter].pairs->erase(column_vector[counter].pairs->begin() +
    choice + column_vector[counter].start);
  column_vector[counter].end--;
  //if (column_vector[counter].start != column_vector[counter].end)
    //column_vector[counter].end--;
  
  if (column_vector[0].top != column_vector[0].pairs->at(column_vector[0].start).sim) {
    column_vector[0].top = column_vector[0].pairs->at(column_vector[0].start).sim;
    for (uint i = 0; i < init_ptr; ++i) {
      column_vector[i].end = bin_search(*(column_vector[i].pairs),
        column_vector[0].top - delta,
        column_vector[i].getStart(),
        column_vector[i].getEnd());
    }
  }
  // if (column_vector[counter].top != 
  //   column_vector[counter].pairs->at(column_vector[counter].start).sim) {
  //   column_vector[counter].top = column_vector[counter].pairs->at(column_vector[counter].start).sim;
  //   column_vector[counter].end = bin_search(*(column_vector[counter].pairs),
  //       column_vector[0].top - delta,
  //       column_vector[i].getStart(),
  //       column_vector[i].getEnd());
  // }
  
  } while (right_exclude.find(to_return.second) != right_exclude.end());

  //while (right_exclude.find(column_vector[result.first].pairs->at(result.second).node) != 
  //right_exclude.end());


  column_vector.erase(column_vector.begin() + left_index);
  init_ptr--;
  for (uint i = 0; i < column_vector.size(); ++i) {
    if (left_exclude.find(column_vector[i].left_node) != left_exclude.end()) {
      if (i < init_ptr)
        init_ptr--;
      column_vector.erase(column_vector.begin() + i);
    }
  }

  //std::sort(column_vector.rbegin(),column_vector.rend());
  //this->debug();
  //exit(0);
  //return std::make_pair(0,0);
  return to_return;
}


std::pair<ushort,ushort> SeedMatrix::pick_random(uint choices){
  std::uniform_int_distribution<int> uni(0, choices - 1);
  std::random_device rd;
  std::mt19937 gen(rd());

  uint chose = uni(gen);
  std::cout << "pick random " << choices << " and " << chose << std::endl;
  uint counter = 0;
  for (uint i = 0; i < column_vector.size(); i++) {
    if ((column_vector[i].end - column_vector[i].start) > (chose - counter)) 
      return std::make_pair(i, chose-counter); 
    else
      counter += (column_vector[i].end - column_vector[i].start);
  }
  return std::make_pair(0,0);
}

uint SeedMatrix::bin_search(vector<node_t> & vec, float val, uint i, uint j){
  std::cout << "bin search " << val << " (" 
	    << i << ", " << j << ") " << std::endl;
  if(i == j || i == j - 1){
    return j;
  }
  /*
  uint m;
  uint k;
  m = i;
  k = j;

  std::cout << "i + j:" << (i + j) << std::endl;
  std::cout << "m + k:" << (m + k) << std::endl;
  m += k;
  
  m /= 2;
  */

  uint m = (i + j) / 2;
  std::cout << "m " << m << " AND " << vec.size() << std::endl;
  std::cout << "sim[m] " << vec[m].sim << std::endl;
  //std::cout << "vector size " << vec.size() << std::endl;
  //std::cout << "before if statements" << std::endl;
  if(vec[m].sim >= val && vec[m+1].sim < val){ 
  //&& ((m + 1) < vec.size()) && vec[m+1].sim > val){
    std::cout << "found " << std::endl;
    return m;
  }else if(vec[m].sim >= val){
    std::cout << "right" << std::endl;
    return bin_search(vec, val, m, j);
  }else{ //vec[m].sim < val
    std::cout << "left" << std::endl;
    return bin_search(vec, val, i, m);
  }
  return 0;
}

vector<node_t> * SeedMatrix::create_node(ushort node_num){
  std::cout << "create node " << node_num << std::endl;
  vector<float> & sim_row = sims[node_num];
  std::cout << "row size " << sim_row.size() << std::endl;
  vector<node_t> * out = new vector<node_t>();
  out->reserve(sim_row.size());
  std::cout << "out vector created" << std::endl;
  for(ushort j = 0; j < sim_row.size(); ++j){
    node_t curr = {sim_row[j], j};
    //std::cout << "node " << j << curr.sim << std::endl;
    out->push_back(curr);
  }
  std::sort(out->rbegin(), out->rend());
  /*
  for(auto it = out->begin(); it != out->end(); ++it){
    std::cout << "(" << it->sim << ", " << it->node << ")" << std::endl;
  }
  */
  std::cout << "create done" << std::endl;
  return out;
}

void SeedMatrix::debug(){
  for(auto it = column_vector.begin(); it != column_vector.end(); ++it){
    seed_t entry = *it;
    std::cout << "( " << entry.left_node << ", " 
	      << ((entry.pairs != NULL) ? "exist" : "Null") << ", "
	      << entry.top << ", " 
	      << entry.start << ", " 
	      << entry.end << ") " << std::endl;
  }
}
