#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <iostream> //purely for testing

template<typename ValueType>
class Trie
{
public:
	Trie();
	~Trie();
	void reset();
	void insert(const std::string& key, const ValueType& value);
	std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

	// C++11 syntax for preventing copying and assignment
	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;
private:
	struct Node
	{
		std::vector<ValueType> vec;
		Node* child[128] = { nullptr };
		std::vector<Node*> ptrVec;
	};
	Node* m_root;
	void cleanUp(Node* current);
	void lookup(std::string key, std::vector<ValueType>& myVec, Node* current, bool offByOneAlreadyUsed, int pos) const;
};


template<typename ValueType>
void Trie<ValueType>::cleanUp(Node* current)
{
	if (current == nullptr)
		return;
	
	typename std::vector<Node*>::iterator it;
	it = current->ptrVec.begin();
	while (it != current->ptrVec.end()) //iterate through vector and recursively call cleanUp
	{
		Node* x = (*it);
		cleanUp(x);
		it++;
	}
	
	delete current;

}

template<typename ValueType>
Trie<ValueType>::Trie()
{
	m_root = new Node;
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
	cleanUp(m_root);
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
	cleanUp(m_root);
	m_root = new Node;
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
	Node* current = m_root; //start current pointer at top of trie

	for (int i = 0; i != key.size(); i++)
	{
		int s = key[i]; //get digit conversion of character in string
		if (current->child[s] == nullptr) //if the node doesn't exist create a new one
		{
			current->child[s] = new Node;
			current->ptrVec.push_back(current->child[s]);
		}
		current = current->child[s];
	}

	current->vec.push_back(value); //place value into current Node's vector
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
	std::vector<ValueType> myVec; //the vector where everything will go 

	Node* current = m_root;
	int s = key[0]; //we know the first character must match

	if (current->child[s] == nullptr)
		return myVec;
	else
	{
		current = current->child[s];
		if(exactMatchOnly)
			lookup(key, myVec, current, true, 1); //in this case the off by one has already been used(since we need an exact match)
		else
			lookup(key, myVec, current, false, 1);
	}
	

	return myVec;
}

template<typename ValueType>
void Trie<ValueType>::lookup(std::string key, std::vector<ValueType>& myVec, Node* current, bool offByOneAlreadyUsed, int pos) const
{

	if (current == nullptr)
		return;

	if (pos == key.size()) //we have now reached the node with the data we want
	{
		typename std::vector<ValueType>::iterator it;
		it = current->vec.begin();
		while (it != current->vec.end())
		{
			myVec.push_back(*it);
			it++;
		}
		return;
	}

	int s = key[pos];
	pos++;
	lookup(key, myVec, current->child[s], offByOneAlreadyUsed, pos); //call lookUp recursively on the correct path
	if (!offByOneAlreadyUsed)
	{
		typename std::vector<Node*>::iterator it;
		it = current->ptrVec.begin();
		while (it != current->ptrVec.end())
		{
			Node* x = (*it);
			if (x != current->child[s])
			{
				lookup(key, myVec, x, true, pos); //call lookUp recursively on all other paths for potential snips
			}
			it++;
		}
	}
}
#endif // TRIE_INCLUDED
