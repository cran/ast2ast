/*
R package etr
Copyright (C) 2021 Konrad Krämer

This file is part of R package etr


etr is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with etr
If not see: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html#SEC4
*/

#ifndef CONCATENATE_H
#define CONCATENATE_H

#include "vec.hpp"

namespace etr {

template<class F, class...Args>
inline F for_each_arg(F f, Args&&...args) {
  (f(std::forward<Args>(args)),...);
  return f;
}

template<typename...Args>
inline VEC<double> coca(Args&&...args) {
  int size = 0;
  for_each_arg([&](auto arg){
    if constexpr (std::is_same<decltype(arg), int>::value) {
      size++;
    } else if(std::is_same<decltype(arg), double>::value) {
      size++;
    } else {
      if constexpr (is_vec<decltype(arg)>::value){
        size += arg.size();
      } 
    }
  }, args...);
  
  VEC<double> ret(size, 0.0);
  int index = 0;
  
  for_each_arg([&](auto arg){
    if constexpr (std::is_same<decltype(arg), int>::value) {
      ret[index] = static_cast<double>(arg);
      index++;
    } else if(std::is_same<decltype(arg), double>::value) {
      ret[index] = arg;
      index++;
    } else {
      if constexpr (is_vec<decltype(arg)>::value){
        for(int i = 0; i < arg.size(); i++) {
          ret[index + i] = arg[i];
        }
        index += arg.size();  
      }
    }
  }, args...);
  
  return ret;
}



/*

 inline VEC<double> combine(VEC<double>&& a, VEC<double>&& b) {
 //a.ismatrix = false;
 //b.ismatrix = false;
 VEC<double> ret(a.size() + b.size());
 for(int i = 0; i < ret.size(); i++) {
 if(i < a.size()) {
 ret[i] = a[i];
 } else if(i >= a.size()) {
 ret[i] = b[i - (a.size())];
 }
 }
 ret.ismatrix = false;
 return ret;
 }
 
 inline VEC<double> combine(VEC<double>& a, VEC<double>& b) {
 //a.ismatrix = false;
 //b.ismatrix = false;
 VEC<double> ret(a.size() + b.size());
 for(int i = 0; i < ret.size(); i++) {
 if(i < a.size()) {
 ret[i] = a[i];
 } else if(i >= a.size()) {
 ret[i] = b[i - (a.size())];
 }
 }
 ret.ismatrix = false;
 return ret;
 }
 
 inline VEC<double> combine(VEC<double>& a, double b) {
 VEC<double> ret(a.size() + 1);
 //a.ismatrix = false;
 for(int i = 0; i < ret.size(); i++) {
 if(i < a.size()) {
 ret[i] = a[i];
 } else if(i >= a.size()) {
 ret[i] = b;
 }
 }
 ret.ismatrix = false;
 return ret;
 }
 
 inline VEC<double> combine(double a, VEC<double>& b) {
 VEC<double> ret(b.size() + 1);
 //b.ismatrix = false;
 for(int i = 0; i < ret.size(); i++) {
 if(i == 0) {
 ret[i] = a;
 } else {
 ret[i] = b[i -1];
 }
 }
 ret.ismatrix = false;
 return ret;
 }
 
 
template <typename ... Ts>
inline VEC<double> coca (Ts && ... multi_inputs) {
    VEC<double> ret;
    ret.ismatrix = false;
    int i = 0;

    ([&] (auto & input)
    {
      if(i == 0) {
        ret = input;
      } else {
        ret = combine(ret, input);
      }

      i++;
    } (multi_inputs), ...);

    return ret;
}




template <typename ... Ts>
inline VEC<double> coca (Ts & ... multi_inputs) {
    VEC<double> ret;
    ret.ismatrix = false;
    int i = 0;

    ([&] (auto & input)
    {
      if(i == 0) {
        ret = input;
      } else {
        ret = combine(ret, input);
      }

      i++;
    } (multi_inputs), ...);

    return ret;
}
*/

}

#endif
