#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
	GenomeImpl(const string& nm, const string& sequence);
	static bool load(istream& genomeSource, vector<Genome>& genomes);
	int length() const;
	string name() const;
	bool extract(int position, int length, string& fragment) const;
private:
	string m_name;
	string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
	// This compiles, but may not be correct
	m_name = nm;
	m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
	
	string s; //our getline string
	string sequence; 
	string name;
	int totalGenomesExtracted = 0;
	while (getline(genomeSource, s))
	{
		if (s.size() == 0) //blank line check
			return false;

		if(s[s.size() - 1] == '\n' || s[s.size() - 1] == '\r') //if we have a newline character 
			s = s.substr(0, s.size() - 1); //create a substring so that newline character is excluded 

		if (s.size() == 0) //blank line check
			return false;

		if (s[0] == '>')
		{
			if (sequence != "" || name != "") //if we're not on our first genome
			{
				if (name == "" || sequence == "") //if a name or sequence was never read in
					return false;
				else
				{
					Genome x(name, sequence); //insert the genome
					genomes.push_back(x);
					sequence = ""; //reset the strings
					name = "";
					totalGenomesExtracted++;
				}
			}

			if (s.size() <= 1) // no characters after > sign
				return false;
			name = s.substr(1, s.size() - 1);
			continue;
		}

		if (s.size() > 80) //more that 80 characters in sequence line
			return false;
		for (int i = 0; i != s.size(); i++)
		{
			switch (s[i])
			{
			case 'a':
			case 'A': sequence += 'A'; break;
			case 't':
			case 'T': sequence += 'T'; break;
			case 'g':
			case 'G': sequence += 'G'; break;
			case 'c':
			case 'C': sequence += 'C'; break;
			case 'n':
			case 'N': sequence += 'N'; break;
			default: return false;
			}
		}
	}

	if (name == "" || sequence == "") //if a name or sequence was never read in for the final genome
		return false;
	else
	{
		Genome x(name, sequence); //insert the genome
		genomes.push_back(x);

		sequence = ""; //reset the strings
		name = "";
		totalGenomesExtracted++;
	}
	return totalGenomesExtracted;
}

int GenomeImpl::length() const
{
	return m_sequence.size();
}

string GenomeImpl::name() const
{
	return m_name;  
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	if (position + length > m_sequence.size())
		return false;
	else
	{
		fragment = m_sequence.substr(position, length);
	}
	return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
	m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
	delete m_impl;
}

Genome::Genome(const Genome& other)
{
	m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
	GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
	delete m_impl;
	m_impl = newImpl;
	return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
	return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
	return m_impl->length();
}

string Genome::name() const
{
	return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
	return m_impl->extract(position, length, fragment);
}
