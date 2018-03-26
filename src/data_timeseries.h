/**
 * @file data_timeseries.h
 * @brief Keeps a timeseries and yields statistic -- average, stddev, min, max.
 * @author Martin Becker <becker@rcs.ei.tum.de>
 * @date 18.04.2014
 
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
 */

#ifndef DATA_TIMESERIES_H
#define DATA_TIMESERIES_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>
#include <algorithm>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <math.h>
#include <iomanip>
#include "data_timed.h"
#include "time_fun.h"


/**
 * @brief a timeseries is a data item whereas each point has a time annotation.
 * This is a template class specific to the data points.
 *
 * TODO: move common stuff to DataTimed
 */
template <typename T>
class DataTimeseries : public DataTimed {
private:
    typedef std::pair<double,T> datapair; ///< one item in the timeline is this

public:
    struct TimedSample{
        double time;
        T data;
        bool operator<(const TimedSample& r) const { return time < r.time; }
        operator datapair() const { return {time,data};}
    };

    /**
     * @brief Statistics
     * @param keepitems if true then individual items are stored.
     */
    DataTimeseries(std::string name, bool keepitems=true) : DataTimed(name), _keepitems(keepitems) {
        _defaults();
    }

    void _defaults() {
        _min_valid = false;
        _max_valid = false;
        _sqsum = _sum = 0.;
        _n = 0;
        _min_t = INFINITY;
        _max_t = -INFINITY;
        _valid = false;
    }

    /**
     * @brief create a new averaged (moving window) dataseries
     * based on the current one.
     * @param windowlen_sec the window length in seconds
     * TODO: optimize
     */
    void moving_average(DataTimeseries<T> & other, float windowlen_sec) const {
        other.clear();

        for (unsigned int k=0; k< this->_elems.size(); ++k) {
            const double t = this->_elems[k].time;
            const double tmin = t-windowlen_sec;
            const double tmax = t+windowlen_sec;
            unsigned n = 0;
            T sum = 0.;
            // left half
            for (unsigned left = k; left > 0; left--) {
                if (this->_elems[left].time < tmin) break;
                n++;
                sum+=this->_elems[left].data;
            }
            // right half
            for (unsigned right = k+1; right < this->_elems.size(); ++right) {
                if (this->_elems[right].time > tmax) break;
                n++;
                sum+=this->_elems[right].data;
            }
            // average window
            if (n > 0) {
                sum /= n;
            }
            other.add_elem(sum, t);
        }
        assert (other._elems.size() == this->_elems.size()); // otherwise something above is wrong
    }

    // implements Data::clear()
    void clear() {
        _defaults();
        _elems.clear();
        _time_epoch_datastart_usec = 0;
    }

    // implements Data::get_typename()
    std::string get_typename(void) const {
        std::string str = "data_timeseries<";
        str += DEBUGTYPE(_min);
        str += ">";
        return str;
    }

    void add_elem(T const &dataelem, double datatime = NAN) {
        _sum += dataelem;
        _sqsum += pow(dataelem,2);
        if (_keepitems) {
            _elems.push_back({datatime,dataelem});
        }
        if (_min_valid) {
            if (dataelem < _min) _min = dataelem;
            if (datatime < _min_t) _min_t = datatime;
        } else {
            _min = dataelem;
            _min_t = datatime;
            _min_valid = true;
        }
        if (_max_valid) {
            if (dataelem > _max) _max = dataelem;
            if (datatime > _max_t) _max_t = datatime;
        } else {
            _max = dataelem;
            _max_t = datatime;
            _max_valid = true;
        }
        _n++;
        _valid = true;
    }

    double get_stddev() const {        
        return sqrt(_sqsum/pow(_sum,2));
    }

    double get_average() const {
        if (_n > 0) {
            return _sum/((double)_n);
        }
        return NAN;
    }

    datapair get_first() const {
        if (_keepitems) {
            if (_n > 0) {
                return _elems.front();
            }
        }
        return datapair(NAN, 0);
    }

    datapair get_last() const {
        if (_keepitems) {
            if (_n > 0) {
                return _elems.back();
            }
        }
        return datapair(NAN, T{});;
    }

    /**
     * @brief _get_index_of_time
     * @param timeinstant time to search for
     * @param idx returns the index, which is >= time
     * @return if true, then idx is valid, otherwise not found
     */
    bool _get_index_of_time(double timeinstant, unsigned int & idx_before, unsigned int & idx_after) const {
        if (timeinstant > _max_t || timeinstant < _min_t) return false; // extrapolation not supported

        // find item which is >= timeinstant
        double t_pre=0.;
        bool first = true;
        for(unsigned int k=0; k< _elems.size(); ++k) {
            const double t = _elems[k].time;
            // assumption: vector is ordered by time
            if (first) {
                t_pre = t;
                first = false;
            }
            //assert(t >= t_pre);
            if (t < t_pre) {
                std::cerr << "time is non-monotonic in data " << this->get_fullname(this) << std::endl;
                return false;
            }
            // --
            if (t == timeinstant) {
                idx_before = k;
                idx_after = k;
                return true;
            } else if (t > timeinstant) {
                // lin. interpolate
                if (k > 0) {
                    idx_before = k - 1;
                    idx_after = k;
                    return true;
                }
            }
            t_pre = t;
        }
        return false;
    }

    /**
     * @brief gets the value of the timeseries that belongs to time instant 'timeinstant'. If there is no data point, it interpolates.
     * @param timeinstant time at which the value shall be retrieved, internal relative time!!
     * @param val value will be written to here
     * @return true if value could be calculated, else false
     */
    bool get_data_at_time(double timeinstant, T &val) const {
        if (timeinstant > _max_t || timeinstant < _min_t) return false; // extrapolation not supported

        unsigned int idx_before = UINT_MAX;
        unsigned int idx_after = 0;
        bool ret = _get_index_of_time(timeinstant, idx_before, idx_after);
        if (!ret) return false;

        assert(idx_before <= idx_after);
        if (idx_before == idx_after) {
            // hit a sample
            val = _elems[idx_before].data;
            return true;
        } else {
            // betweem samples; need interpolation
            TimedSample pre = _elems[idx_before];
            TimedSample post = _elems[idx_after];
            const double m = (post.data-pre.data)/(post.time-pre.time);
            val = pre.data + (timeinstant - pre.time)*m;
            return true;
        }
    }
    TimedSample get_data_at_time(double timeinstant) const {
        TimedSample sample{};
        sample.time = timeinstant;
        auto success =  get_data_at_time(timeinstant, sample.data);
        if (!success) {
            std::stringstream ss;
            ss << "Could not find or interpolate data for time instant" << timeinstant
               << " in timed series \"" << get_fullname <<"\"";
            throw std::out_of_range(ss.str());
        }
        return sample;
    }

    /**
     * @brief get_min_time
     * @return minimum TIME in data series
     */
    double get_min_time() const {
        return _min_t;
    }

    /**
     * @return maximum TIME in data series
     */
    double get_max_time() const {
        return _max_t;
    }

    /**
     * @brief get_min
     * @return  mininum VALUE in data series
     */
    T get_min() const {
        return _min;
    }



    bool get_stats_timewindow(double tmin, double tmax, data_stats & s) const {
        if (tmax < tmin) return false;

        if (tmin < _min_t) tmin = _min_t;
        if (tmax > _max_t) tmax = _max_t;

        // find first at tmin
        unsigned int idx_min_pre, idx_min_post;
        bool ret = _get_index_of_time(tmin, idx_min_pre, idx_min_post);
        if (!ret) return false;

        // find last at tmax
        unsigned int idx_max_pre, idx_max_post;
        ret = _get_index_of_time(tmax, idx_max_pre, idx_max_post);
        if (!ret) return false;

        // we have the limits
        s.n_samples = idx_max_pre - idx_min_post + 1;

        double sum = 0.;
        double sumsq = 0.;
        double first = true;
        unsigned int n_samples_int = 0; // internally we may use more samples. The returned struct says "samples in between", but we may interpolate at beginning and end.

        // do we need to interpolate at the beginning?
        if (idx_min_pre != idx_min_post) {
            // yes
            T val;
            if (get_data_at_time(tmin, val)) {
                s.min = (double)val;
                s.max = (double)val;
                first = false;
                sum += val;
                sumsq += (val*val);
                n_samples_int++;
            }
        }

        // run over the samples now
        for (unsigned int k = idx_min_post; k <= idx_max_pre; k++) {
            const T val = _elems[k].data;
            if (first) {
                s.min = (double)val;
                s.max = (double)val;
                first = false;
            } else {
                s.min = (val < s.min) ? val : s.min;
                s.max = (val > s.max) ? val : s.max;
            }
            sum += val;
            sumsq += (val*val);
            n_samples_int++;
        }

        // do we need to interpolate at the end?
        if (idx_max_pre != idx_max_post) {
            // yes
            T val;
            if (get_data_at_time(tmax, val)) {
                if (first) {
                    s.min = val;
                    s.max = val;
                    first = false;
                } else {
                    s.min = (val < s.min) ? val : s.min;
                    s.max = (val > s.max) ? val : s.max;
                }
                sum += val;
                sumsq += (val*val);
                n_samples_int++;
            }
        }

        // finally: complete stats
        sum /= n_samples_int;
        sumsq /= n_samples_int;
        s.avg = sum;
        s.stddev = sqrt(sumsq - sum*sum);
        s.t_min = tmin;
        s.t_max = tmax;
        const double dt = tmax - tmin;
        s.freq = dt != 0.0 ? n_samples_int / (dt) : 0.0;
        return true;
    }

    /**
     * @brief get_max
     * @return maximum VALUE of data series
     */
    T get_max() const {
        return _max;
    }

    unsigned long get_epoch_dataend() const {
        if (_elems.empty()) { return get_epoch_datastart(); }
        return ((unsigned long) (_elems.back().time*1E6)) + get_epoch_datastart();
    }

    std::vector<double> get_time() const {
        std::vector<double> ret(_elems.size());
        for (std::size_t i =0; i< _elems.size();++i){
            ret[i] = _elems[i].time;
        }
        return ret;
    }

    std::vector<T> get_data() const {
        std::vector<T> ret(_elems.size());
        for (std::size_t i =0; i< _elems.size();++i){
            ret[i] = _elems[i].data;
        }
        return ret;
    }

    const std::vector<TimedSample>& get_timed_data() const {
        return _elems;
    }

    /**
     * @brief fetch data at given index
     * @return true if fetched, else data is invalid
     */
    bool get_data(unsigned int index, double &tval, T &dval) const {
        if (index > _n) return false;
        dval = _elems[index].data;
        tval = _elems[index].time;
        return true;
    }

    TimedSample get_data(unsigned int index) const {
        if (index > _n) {
            std::stringstream ss;
            ss << "Tried to access " << index <<" th element "
               << "in timed series \"" << get_fullname <<"\n ,"
               << "which only has " << _n << " elements.";
            throw std::out_of_range(ss.str());
        }
        return _elems[index];
   }

    /**
     * @brief how many elements are stored?
     * @return
     */
    unsigned int size() const {
        return _n;
    }

    /**
     * @return rate of data (in Hz)
     */
    float get_rate() const {
        const double timespan = get_max_time() - get_min_time();
        float freq = 0.0;
        if (timespan > 0.0) freq = size() / timespan;
        return freq;
    }

    std::string describe_myself() const {
        std::stringstream ss;
        ss<< "Data: TIMESERIES" << std::endl <<
             "type: " << get_typename() << ", " << _verbose_data_class() << std::endl <<
             "name: " << Data::get_fullname(this) << std::endl <<
             "units:" << _units << std::endl <<             
             "#data points: " << _n << std::endl <<
             "data rate: ~" << get_rate() << "Hz" << std::endl <<
             "time: "<< get_min_time() << " .. " << get_max_time() << std::endl <<
             "min: " << get_min() << " " << _units << std::endl <<
             "max: " << get_max() << " " << _units << std::endl <<
             "avg: " << get_average() << " " << _units  << std::endl <<
             "stddev: " << get_stddev() << " " << _units  << std::endl <<
             "time bad: " << (has_bad_timestamps() ? "true" : "false") << std::endl;

        return ss.str();
    }

    /**
     * @brief distributes all samples over the time span equidistantly.
     */
    void make_periodic(void) {
        double rate = get_rate();
        if (rate <= 0.0) return;

        double dt = 1.0 / get_rate();
        double t0 = get_min_time();
        for (unsigned int k=0; k<_elems.size(); ++k) {
            _elems[k].time = t0 + dt*k;
        }

        _bad_timestamps = false;
        _class = DATA_DERIVED;
    }

    // implements Data::export_csv()
    bool export_csv(const std::string & filename, const std::string & sep = std::string(",")) const {
        std::ofstream fout(filename.c_str());

        if (!fout.is_open()) return false;

        fout << "#time, " << _name << "[" << _units << "]" << std::endl;
        fout << std::setprecision(9);
        for (unsigned int k=0; k<_n; k++) {
            // write
            fout << _elems[k].time << sep << _elems[k].data << std::endl;
        }

        fout.close();
        return true;
    }

    /**
     * @brief implements Data::Clone()
     * @return covariant of Data*
     */
    DataTimeseries* Clone() const {
        return new DataTimeseries(*this);
    }

    // implements Data::merge_in()
    bool merge_in(const Data * const other) {
        const DataTimeseries*const src = dynamic_cast<const DataTimeseries*const>(other);
        if (!src) return false;
        if (!src->_valid) return false;

        const double tmin_src = src->_elems.front().time;
        const double tmax_src = src->_elems.back() .time;
        const double tmin_me = _elems.front().time;
        const double tmax_me = _elems.back().time ;

        /*****************
         *  MERGING IN
         *****************/

        if (tmax_src < tmin_me) {
            _elems.insert(_elems.begin(), src->_elems.begin(), src->_elems.end()); ///< prepend data
        } else if (tmin_src > tmax_me){
            _elems.insert(_elems.end(), src->_elems.begin(), src->_elems.end()); ///< append data
        } else {
            // Merge: data is overlapping...we have to sort-in every single data item

            const auto& own = _elems;
            const auto& others = src->_elems;

            //c++11:
            assert(std::is_sorted(own.begin(),own.end()));
            assert(std::is_sorted(others.begin(),others.end()));

            //merge into new buffer
            const size_t total_size = own.size() + others.size();
            std::vector<TimedSample> merged(total_size);
            std::merge(
                        own.begin(),
                        own.end(),
                        others.begin(),
                        others.end(),
                        merged.begin()
                        );
            _elems = std::move(merged);
        }

        // correct the other class members
        _sum+=src->_sum;
        _sqsum+=src->_sqsum;
        _max= std::max(_max,src->_max);
        _min= std::min(_min,src->_min);
        _min_t = _elems.front().time;
        _max_t = _elems.back().time;
        _n+= src->_elems.size();

        return true;
    }

private:
    std::size_t    _n;
    bool            _keepitems;

    std::vector<TimedSample> _elems;

    double          _sum;
    double          _sqsum;
    T               _max;
    T               _min;
    double          _max_t;
    double          _min_t;
    bool            _max_valid;
    bool            _min_valid;
};

#endif // DATA_TIMESERIES_H
