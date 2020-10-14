#include "provided.h"
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Trie.h"
using namespace std;

class GenomeMatcherImpl
{
public:
	GenomeMatcherImpl(int minSearchLength);
	~GenomeMatcherImpl();
	void addGenome(const Genome& genome);
	int minimumSearchLength() const;
	bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
	bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	struct TrieMatch
	{
		Genome* ptr;
		int position;
	};

	int m_minSearchLength;
	Trie<TrieMatch> myTrie;
	vector<Genome*> genomeVec;
	void addToMatchVec(string name, int pos, int size, int minimumLength, vector<DNAMatch>& matches) const;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
	m_minSearchLength = minSearchLength;
}

GenomeMatcherImpl::~GenomeMatcherImpl()
{
	vector<Genome*>::iterator it = genomeVec.begin();
	while (it != genomeVec.end())
	{
		delete (*it);
		it++;
	}
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	if (m_minSearchLength > genome.length())
		return;

	Genome* x = new Genome(genome);
	genomeVec.push_back(x);

	for (int i = 0; i <= genome.length() - m_minSearchLength; i++)
	{
		TrieMatch y;
		y.position = i;
		y.ptr = x;
		string s;
		if(genome.extract(i, m_minSearchLength, s))
			myTrie.insert(s, y);
	}
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return m_minSearchLength;
}

void GenomeMatcherImpl::addToMatchVec(string name, int pos, int size, int minimumLength, vector<DNAMatch>& matches) const
{
	if (size >= minimumLength)
	{
		vector<DNAMatch>::iterator it = matches.begin();
		while (it != matches.end()) 
		{
			if ((*it).genomeName == name) //if the genomeName already exists
			{
				if ((*it).length < size) //keep the longer sequence in matches vector
				{
					(*it).length = size;
					(*it).position = pos;
					return;
				}
				else if ((*it).length > size)
					return;
				else if ((*it).length == size) //if the sequences match in length
				{
					if ((*it).position < pos) //keep the earlier position in matches vector
						return;
					else
					{
						(*it).position = pos;
						return;
					}	
				}
			}
			it++;
		}

		DNAMatch x;
		x.genomeName = name;
		x.position = pos;
		x.length = size;
		matches.push_back(x);
	}
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if (fragment.size() < minimumLength || minimumLength < m_minSearchLength)
		return false;

	vector<DNAMatch>::iterator it1 = matches.begin(); //empty out the matches vector because we will now be using it
	while (it1 != matches.end())
		it1 = matches.erase(it1);
	
	string s = fragment.substr(0, m_minSearchLength);
	vector<TrieMatch> finds = myTrie.find(s, exactMatchOnly);
	
	
	vector<TrieMatch>::iterator it = finds.begin();
	while (it != finds.end())
	{
		int offbyOneCounter = 0;
		int pos = (*it).position;
		for (int i = 0; i != fragment.size(); i++) 
		{
			string potential;
			if (!((*it).ptr->extract(pos, i + 1, potential)))
				break;

			bool reachedEnd = false;

			if ((*it).ptr->length() == pos + i + 1)
				reachedEnd = true;

			if (potential == fragment) //if we have a perfect match
			{
				addToMatchVec((*it).ptr->name(), pos, potential.size(), minimumLength, matches);
				break;
			}

			if (exactMatchOnly)
			{
				if (potential != fragment.substr(0, i + 1) || reachedEnd) //if we reach a point where the strings mismatch or the end
				{
					if(reachedEnd)
						addToMatchVec((*it).ptr->name(), pos, potential.size(), minimumLength, matches);
					else
						addToMatchVec((*it).ptr->name(), pos, i, minimumLength, matches);	
					break;
				}
			}
			else
			{
				if (potential != fragment.substr(0, i + 1)) //if we reach a point where the strings mismatch 
				{
					if (potential.at(potential.size() - 1) != fragment.at(i))
						offbyOneCounter++;
					if (reachedEnd && offbyOneCounter != 2)
					{
						addToMatchVec((*it).ptr->name(), pos, potential.size(), minimumLength, matches);
						break;
					}
					if (offbyOneCounter == 2) //if we are off by 2
					{
						addToMatchVec((*it).ptr->name(), pos, i, minimumLength, matches);
						break;
					}

					if (i + 1 == fragment.size())
					{
						addToMatchVec((*it).ptr->name(), pos, fragment.size(), minimumLength, matches);
						break;
					}
				}
			}

		}
		it++;
	}
	return !matches.empty();
}

bool querySort(const GenomeMatch& a, const GenomeMatch& b)
{
	if (a.percentMatch > b.percentMatch) //sort based on greater percentMatch
		return true;
	else if (a.percentMatch < b.percentMatch)
		return false;

	return !(a.genomeName > b.genomeName); //sort based on alphabetical order of names
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	if (fragmentMatchLength < m_minSearchLength)
		return false;


	int s = query.length() / fragmentMatchLength;
	map<string, int> name2count; //a map that tracks genome names and the number of times they have appeared

	for (int i = 0; i != s; i++)
	{
		vector<DNAMatch> match;
		string extract;
		query.extract(i*fragmentMatchLength, fragmentMatchLength, extract); //extract string from query
		findGenomesWithThisDNA(extract, fragmentMatchLength, exactMatchOnly, match); //find matches and store in match vector

		vector<DNAMatch>::iterator it = match.begin(); //iterate through match vector
		while (it != match.end())
		{
			map<string, int>::iterator it1 = name2count.find((*it).genomeName); //search for the genomeName in map
			if (it1 == name2count.end()) //if the genomeName is not found in our map, add new element to map
			{
				name2count[(*it).genomeName] = 1;
			}
			else //if the genomeName is found in our map, increment the count
			{
				(*it1).second++;
			}
			it++;
		}
	}

	vector<GenomeMatch>::iterator resultsIterator = results.begin(); //empty out results vector becuase at this point we are using it 
	while (resultsIterator != results.end())
		resultsIterator = results.erase(resultsIterator);

	map<string, int>::iterator mapIterator = name2count.begin(); //iterate through map
	while (mapIterator != name2count.end())
	{
		double percentage = ((*mapIterator).second * 100.0) / s; //calculate percent match 
		if (percentage >= matchPercentThreshold) //add to results if percentage meets threshold requirement
		{
			GenomeMatch x;
			x.genomeName = (*mapIterator).first;
			x.percentMatch = percentage;
			results.push_back(x);
		}
		mapIterator++;
	}

	sort(results.begin(), results.end(), &querySort); //sort results vector using custom sorting function 

	return !results.empty();
}



//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
	m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
	delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
	m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
	return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}