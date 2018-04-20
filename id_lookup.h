//
// Created by Martin Rusilowicz on 17/08/2017.
//

#ifndef COINFINDER_IDLOOKUP_H
#define COINFINDER_IDLOOKUP_H

#include <map>
#include <string>
#include <unordered_set>
#include <set>



/**
 * Forms a simple lookup-by-name dictionary that creates elements if they don't already exist.
 * @tparam T    Type of element. 
 */
template <typename T>
class id_lookup
{
    private:
        std::map<std::string, T*>* _id_lookup;
        
    public:
        id_lookup();
        ~id_lookup();
        T& find_id( const std::string& name ) const;
        const std::map<std::string, T*>& get_table() const;
        std::map<std::string, T*>& get_table();
        int size() const;
};


#endif //COINFINDER_IDLOOKUP_H
