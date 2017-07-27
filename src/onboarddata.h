/**
 * @file onboarddata.h
 * @brief This class represents a data item that is logged to the onboard data flash
 * @author Martin Becker <becker@rcs.ei.tum.de>
 * @date 11/12/2014
 
    This file is part of MavLogAnalyzer, Copyright 2014 by Martin Becker.
    
    MavLogAnalyzer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
    
 */

#ifndef ONBOARDDATA_H
#define ONBOARDDATA_H

#include <inttypes.h>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>


class OnboardData
{
    /**********************************
     *  TYPES
     **********************************/
    template<class K, class T>
    struct FlatMap {
        using mapped_type = T;
        using key_type = K;

    private:
        using Entry_t = std::pair<key_type, mapped_type>;
        using Store_t = std::vector<Entry_t>;

    public:
        using iterator       = typename Store_t::iterator;
        using const_iterator = typename Store_t::const_iterator;

        iterator begin()    { return items.begin();}
        iterator end()      { return items.end();}
        const_iterator begin()  const { return items.begin();}
        const_iterator end()    const { return items.end();}

        iterator find(const std::string& s) {
            return std::find_if(items.begin(),items.end(),[&](const Entry_t& e){ return s == e.first; });
        }
        const_iterator find(const std::string& s) const {
            return std::find_if(items.begin(),items.end(),[&](const Entry_t& e){ return s == e.first; });
        }

        mapped_type& operator[](const std::string& s) {
            auto it = find(s);
            if ( it != end() ) {
                return it->second;
            } else {
                items.emplace_back(s,T{});
                return items.back().second;
            }
        }

    private:
        Store_t items;
    };

public:
    typedef FlatMap<std::string, float>         floatdata_t;
    typedef FlatMap<std::string, int64_t>       intdata_t;
    typedef FlatMap<std::string, uint64_t>      uintdata_t;
    typedef FlatMap<std::string, bool>          booldata_t;
    typedef FlatMap<std::string, std::string>   stringdata_t;





    /**********************************
     *  METHODS
     **********************************/
    OnboardData() : _valid(false) {}

    /**
     * @brief only if this returns true, this class carries actual data
     * @return
     */
    bool is_valid(void) const { return _valid; }

    /**
     * @brief returns a string describing the collection of data this class contains
     * @return string
     */
    std::string get_message_origname(void) const { return _msgname_orig; }
    std::string get_message_name(void) const { return _msgname_readable; }

    /**
     * @brief access the data
     * @return
     */
    const floatdata_t & get_floatdata(void) const { return _float_data; }
    const intdata_t & get_intdata(void) const { return _int_data; }
    const uintdata_t & get_uintdata(void) const { return _uint_data; }
    const booldata_t & get_booldata(void) const { return _bool_data; }
    const stringdata_t & get_stringdata(void) const { return _string_data; }


public:
    /**********************************
     *  VARIABLES (that should be private)
     **********************************/
    bool _valid;

    std::string  _msgname_orig;  ///< msg name as in the log, e.g. "CTUN"
    std::string  _msgname_readable; ///< a more readable string, e.g. "Controller Tuning"

    // colums by type. first=col name, second=value
    floatdata_t  _float_data;
    intdata_t    _int_data;
    uintdata_t   _uint_data;
    booldata_t   _bool_data;
    stringdata_t _string_data;

};

#endif // ONBOARDDATA_H
