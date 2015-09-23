#ifndef GRAPH_COMMON
#define GRAPH_COMMON

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>


// http://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using std::cout;
using std::endl;
using std::cerr;

using std::vector;
using std::map;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::ifstream;
using std::ofstream;

/**
 * @brief num_pairs
 * @param x
 * @return
 */
inline size_t num_pairs(size_t x)
{
    return x*(x-1)/2;
}


/**
 * @brief mapvalue_sum
 * @param m
 * @return
 */
template <class T>
T mapvalue_sum(const map<size_t,T> &m)
{
    T sum = 0;
    for ( typename map<size_t,T>::const_iterator i = m.begin(); i!=m.end(); ++i)
        sum+= i->second;
    return sum;
}

/**
 * @brief operator <<
 * @param os
 * @param m
 * @return
 */
template <class T>
inline ostream &operator<<(ostream &os, const map<size_t,T> &m)
{
    for ( typename map<size_t,T>::const_iterator i=m.begin() ; i!= m.end() ; i++)
    {
        os << i->first << "-> " << i->second << endl;
    }
    return os;
}


/**
 * @brief operator <<
 * @param os
 * @param x
 * @return
 */
template <class T>
inline ostream& operator<<(ostream& os, vector<T> x)
{
    os << "[ ";
    for ( typename vector<T>::const_iterator it = x.begin(); it!=x.end(); ++it)
        os << *it << ", ";
    os << "]" << endl;
    return os;
}

/**
 * @brief sort_second
 * @param l
 * @param r
 * @return
 */
template <class T>
inline bool sort_second( const pair<int,T>& l,  const pair<int,T>& r)
{
    return l.second > r.second;
}

#endif // COMMON

