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
        bool operator<(const TimedSample& r){ return time < r.time; }
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

        for (unsigned int k=0; k< this->_elems_time.size(); ++k) {
            const double t = this->_elems_time[k], tmin = t-windowlen_sec, tmax=t+windowlen_sec;
            unsigned n = 0;
            double sum = 0.;
            // left half
            for (unsigned left = k; left > 0; left--) {
                if (this->_elems_time[left] < tmin) break;
                n++;
                sum+=this->_elems_data[left];
            }
            // right half
            for (unsigned right = k+1; right < this->_elems_data.size(); ++right) {
                if (this->_elems_time[right] > tmax) break;
                n++;
                sum+=this->_elems_data[right];
            }
            // average window
            if (n > 0) {
                sum /= n;
            }
            other.add_elem(sum, t);
        }
        assert (other._elems_data.size() == this->_elems_data.size()); // otherwise something above is wrong
    }

    // implements Data::clear()
    void clear() {
        _defaults();
        _elems_data.clear();
        _elems_time.clear();
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
            _elems_data.push_back(dataelem);
            _elems_time.push_back(datatime);
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
                return datapair(_elems_time.front(),_elems_data.front());
            }
        }
        return datapair(NAN, 0);
    }

    datapair get_last() const {
        if (_keepitems) {
            if (_n > 0) {
                return datapair(_elems_time.back(),_elems_data.back());
            }
        }
        return datapair(NAN, 0);;
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
        for(unsigned int k=0; k< _elems_time.size(); ++k) {
            const double t = _elems_time[k];
            // assumption: vector is ordered by time
            if (first) {
                t_pre = t;
                first = false;
            }
            //assert(t >= t_pre);
            if (!t >= t_pre) {
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
            val = _elems_data[idx_before];
            return true;
        } else {
            // betweem samples; need interpolation
            const double val_pre = _elems_data[idx_before];
            const double val_post = _elems_data[idx_after];
            const double t_post = _elems_time[idx_after];
            const double t_pre = _elems_time[idx_before];
            const double m = (val_post-val_pre)/(t_post-t_pre);
            val = val_pre + (timeinstant - t_pre)*m;
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

        // convert the absolute time to time of this series
        const double time_offset = _time_epoch_datastart_usec /1E6;
        tmin -= time_offset;
        tmax -= time_offset;
        // now we got internal time
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
            const T val = _elems_data[k];
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
        if (_elems_time.empty()) { return get_epoch_datastart(); }
        return ((unsigned long) _elems_time.back()*1E6) + get_epoch_datastart();
    }

    const std::vector<double>& get_time() const {
        return _elems_time;
    }

    const std::vector<T>& get_data() const {
        return _elems_data;
    }

    /**
     * @brief fetch data at given index
     * @return true if fetched, else data is invalid
     */
    bool get_data(unsigned int index, double &tval, T &dval) const {
        if (index > _n) return false;
        dval = _elems_data[index];
        tval = _elems_time[index];
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
        return { _elems_time[index], _elems_data[index]};
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
             "name: " << Data::get_fullname(dynamic_cast<const Data*const>(this)) << std::endl <<
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
        for (unsigned int k=0; k<_elems_time.size(); ++k) {
            _elems_time[k] = t0 + dt*k;
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
            T d = _elems_data[k];
            double t = _elems_time[k];
            // write
            fout << t << sep << d << std::endl;
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

        const double tmin_src = src->_elems_time.front() + src->_time_epoch_datastart_usec/1E6;
        const double tmax_src = src->_elems_time.back() + src->_time_epoch_datastart_usec/1E6;
        const double tmin_me = _elems_time.front() + _time_epoch_datastart_usec/1E6;
        const double tmax_me = _elems_time.back() + _time_epoch_datastart_usec/1E6;
        const double dt_sec = (_time_epoch_datastart_usec/1E6 - src->_time_epoch_datastart_usec / 1E6); ///< positive, if my data is more recent

        /*
         * we want no negative time stamps, so adjust all *my* relative times by
         * applying the offset between src and me to *me*
         */
        if (dt_sec > 0.) {
            _time_epoch_datastart_usec = src->_time_epoch_datastart_usec;
            for (std::vector<double>::iterator it = _elems_time.begin(); it != _elems_time.end(); ++it) {
                *it += dt_sec;  ///< correct my relative times
            }
        }

        /*****************
         *  MERGING IN
         *****************/
        _elems_time.reserve(_elems_time.size()+src->_elems_time.size());
        _elems_data.reserve(_elems_data.size()+src->_elems_data.size());
        const bool do_fast_merge = (tmax_src < tmin_me) || (tmin_src > tmax_me); ///< checks for non-overlapping time ranges
        if (do_fast_merge) {
            if (dt_sec > 0.) {
                // PREPEND: my data is later (other earlier). we want no negative time stamps, so adjust all my relative times by adding apply the offset from src
                _elems_time.insert(_elems_time.begin(), src->_elems_time.begin(), src->_elems_time.end()); ///< prepend time
                _elems_data.insert(_elems_data.begin(), src->_elems_data.begin(), src->_elems_data.end()); ///< prepend data
            } else {
                // APPEND: my data is older (other more recent). adjust other data's time relative time stamps by adding the offset to it
                for (std::vector<double>::const_iterator it = src->_elems_time.begin(); it != src->_elems_time.end(); ++it) {
                    _elems_time.push_back(*it - dt_sec); ///< correct other's time stamp and append at the same time
                }
                _elems_data.insert(_elems_data.end(), src->_elems_data.begin(), src->_elems_data.end()); ///< append data
            }
        } else {
            // INSERT: data is overlapping...we have to sort-in every single data item

            //merge time and data arrays into one array (SOA to AOS) for both our data and other data
            std::vector<TimedSample> own(_elems_time.size());
            for (size_t cnt = 0; cnt < _elems_time.size(); cnt++) {
                own[cnt] = TimedSample{_elems_time[cnt], _elems_data[cnt]};
            }

            std::vector<TimedSample> others(src->_elems_time.size());
            for (size_t cnt = 0; cnt < src->_elems_time.size(); cnt++) {
                others[cnt] = TimedSample{src->_elems_time[cnt] - dt_sec, src->_elems_data[cnt]};
            }

            //merge the two series

            //c++11:
            //    assert(std::is_sorted(own.begin(),own.end()),"Other data is not sorted");
            //    assert(std::is_sorted(others.begin(),others.end()),"Other data is not sorted");

            const size_t total_size = _elems_time.size() + src->_elems_time.size();
            std::vector<TimedSample> merged(total_size);
            std::merge(
                        own.begin(),
                        own.end(),
                        others.begin(),
                        others.end(),
                        merged.begin()
                        );

            //split array containing time and data back into separate arrays
            _elems_time.resize(total_size);
            _elems_data.resize(total_size);

            for (std::size_t i = 0; i < total_size; ++i) {
                _elems_time[i] = merged[i].time;
                _elems_data[i] = merged[i].data;
            }
        }

        // correct the other class members
        _sum+=src->_sum;
        _sqsum+=src->_sqsum;
        if (src->_max > _max) _max = src->_max;
        if (src->_min < _min) _min = src->_min;
        _min_t = _elems_time.front();
        _max_t = _elems_time.back();
        _n+= src->_elems_data.size();

        return true;
    }

private:
    unsigned int    _n;
    bool            _keepitems;    

    // FIXME: this storage format is not all that good...pairs would be nicer, but are harder to access
    //std::vector< timeseries_elem >  _elems; // better but plotting would be tedious
    std::vector<T>      _elems_data;    ///< only used if keepitems=true
    std::vector<double> _elems_time;    ///< only used if keepitems=true

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
