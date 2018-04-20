//
// Created by Martin Rusilowicz on 17/08/2017.
//

#include "id_lookup.h"
#include "elements.h"


/**
 * Constructor
 * @tparam T 
 */
template <typename T>
id_lookup<T>::id_lookup()
: _id_lookup(new std::map<std::string, T*>())
{
    // pass
}

/**
 * Destructor
 * All elements are freed
 * @tparam T 
 */
template <typename T>
id_lookup<T>::~id_lookup()
{
    for (auto it : *this->_id_lookup)
    {
        delete it.second;
    }
    
    delete this->_id_lookup;
}

/**
 * Finds the element with the specified `name`, creating it if it doesn't already exist. 
 * @tparam T 
 * @param name 
 * @return 
 */
template <typename T>
T& id_lookup<T>::find_id( const std::string& name ) const
{
    auto it = this->_id_lookup->find(name);
    
    if (it == this->_id_lookup->end())
    {
        T* id = new T(name);
        (*this->_id_lookup)[name] = id;
        return *id;
    }
    else
    {
        return *it->second;
    }
}

/**
 * Obtains the underlying map.
 * @tparam T 
 * @return 
 */
template <typename T>
const std::map<std::string, T*>& id_lookup<T>::get_table() const
{
    return *this->_id_lookup;
}

/**
 * Obtains the underlying map.
 * @tparam T 
 * @return 
 */
template <typename T>
std::map<std::string, T*>& id_lookup<T>::get_table()
{
    return *this->_id_lookup;
}

/**
 * Returns the number of elements in the map.
 * @tparam T 
 * @return 
 */
template <typename T>
int id_lookup<T>::size() const
{
    return static_cast<int>(this->_id_lookup->size());
}

// We only use `id_lookup` on these types, so register them now
template class id_lookup<Alpha>;
template class id_lookup<Beta>;
template class id_lookup<Gamma>;
