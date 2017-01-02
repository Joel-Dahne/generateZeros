/*   File: mult_list.h 
 
     A template version of a linked list, 
     with dynamic memory allocation.

     Usage: List<int>, List<vector>, etc...
 
     Latest edit: Sat Nov 25 2000

     Altered by T. Johnson tue 4 april 19.45


*/



#ifndef MULT_LIST_H
#define MULT_LIST_H 

#include <iostream>
//#include <cstdio>
using std::cout;
using std::cerr;
using std::ostream;
using std::endl;

#undef AT_HOME
#define AT_HOME

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
class List;

template<class TYPE>
class ListObject
{
 public:
  friend class List<TYPE>;
  ListObject  *next;
  TYPE  element;
  ListObject () {};
  ListObject(ListObject<TYPE> *n, TYPE e) : next(n), element(e) {}
};

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
ListObject<TYPE> *copyObject(ListObject<TYPE> *p);


template<class TYPE>
class List 
{
 private:
  ListObject<TYPE> *start;	 
  ListObject<TYPE> *end;	 
  ListObject<TYPE> *current;	 
  ListObject<TYPE> *lastcur;	 
  int  len;		 
  int  maxlen;		 
public:
  List()
    { start = end = current = lastcur = __null ; len = maxlen = 0; }
  ~List();
  List(const List & Li);
  //List(List & li){ printf("List (Copy): only references allowed\n"); exit(1); }
  friend int  Finished    (const List & li) { return (li.current == __null); }
  friend int  IsEmpty     (const List & li) { return (li.start   == __null ); }
  friend int  Length      (const List & li) { return li.len; }
  friend int  MaxLength   (const List & li) { return li.maxlen; }
  friend void ResetLength (      List & li) { li.maxlen = li.len; }
  void operator +=  (const TYPE &);
  void operator *=  (const TYPE &);
  void operator --  ();
  const List<TYPE> & operator = (const List<TYPE> &Li); 
  


#if defined(__linux) && !defined(AT_HOME)
  friend TYPE & First   <TYPE>(List &);
  friend TYPE & Next    <TYPE>(      List &);
  friend TYPE & Last    <TYPE>(const List &);
  friend TYPE & Current <TYPE>(const List &);
  friend TYPE   Pop     <TYPE>(      List &);
  friend void RemoveCurrent <TYPE>(List &);
  friend void Clear <TYPE>(List &);
  friend ostream & operator << <TYPE>(ostream &, const List &);
#endif // __linux

  
#if defined(__sparc) || defined(AT_HOME)
  friend TYPE & First   (List & li)
  {
    li.current = li.start;
    li.lastcur = __null;
    if (li.current == __null) 
      { printf("List (First): empty list"); exit(1); }
    return li.current->element;
  }
  friend TYPE & Next    (List & li)
  {
    if ((li.current == __null) && (li.lastcur == __null)) return First(li);
    if (li.current == __null) return li.lastcur->element;
    li.lastcur = li.current;
    li.current = li.current->next;
    return (li.current == __null) ? li.lastcur->element : li.current->element;
  }
  friend TYPE & Last    (const List & li)
  {
    if (li.end == __null) 
      { printf("List (Last): empty list"); exit(1); }
    return li.end->element;
  }
  friend TYPE & Current (const List & li)
  {
    if (li.current == __null) 
      { printf("List (Current): no element"); exit(1); }
    return li.current->element;
  }
  friend TYPE Pop (List & li)
  {
    TYPE r = First(li);
    RemoveCurrent(li);
    return r;
  }
  friend void RemoveCurrent (List & li)
  {
    ListObject<TYPE> *del_cur = li.current;
    if (li.current == __null) 
      { printf("List (RemoveCurrent): no element"); exit(1); }
    li.current = li.current->next;
    delete del_cur;
    li.len--;
    if (li.lastcur == __null) li.start = li.current;
    else li.lastcur->next = li.current;
    if (li.current == __null) li.end = li.lastcur;
  }
  friend void Clear (List &li)
  {
    if( IsEmpty(li) ) return;
    while( !IsEmpty(li) )
      --li;
  }
  friend ostream & operator << (ostream & o, const List & li)
  {
    int  i = 0;
    if (li.start == __null) o << "*EMPTY*\n";
    else {
      for (ListObject<TYPE> *p = li.start; p != __null; p = p->next) {
	o.width(3);
	  o << (++i) << ": ";
	  o << (p->element) << '\n';
      }
    }
    return o;
  }
#endif // __sparc
};

  




////////////////////////////////////////////////////////////////////////////

template<class TYPE>
List<TYPE>::List(const List & Li)  {
  len=Li.len;
  maxlen=len;
  start = copyObject<TYPE>(Li.start);  
  end = copyObject<TYPE>(Li.end);
  current=start;
  lastcur=__null;
}

template<class TYPE>
List<TYPE>::~List()		 
{
  ListObject<TYPE> *temp;
  while (start != __null) {
    temp = start; start = start->next;
    delete temp;
  }
  start = end = __null;
  len = 0;
}


////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void List<TYPE>::operator += (const TYPE & obj)
{
  ListObject<TYPE> *p = new ListObject<TYPE>;
  p->element = obj;
  p->next = __null;
  if (start == __null) start = p;
  else end->next = p;
  end = p;
  len++;
  if (len > maxlen) maxlen = len;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void List<TYPE>::operator *= (const TYPE & obj)
{
  ListObject<TYPE> *p = new ListObject<TYPE>;
  p->element = obj;
  p->next = start;
  if (start == __null) end = p;
  start = p;
  len++;
  if (len > maxlen) maxlen = len;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void List<TYPE>::operator -- ()
{
  if (start == __null) 
    { printf("List (--): empty list\n"); exit(1); }
  ListObject<TYPE> *p = start->next;
  delete start;
  start = p;
  if (start == __null) end = __null;
  len--;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
ostream & operator << (ostream & o, const List<TYPE> &li)
{
  int  i = 0;
  if (li.start == __null) o << "*EMPTY*\n";
  else {
    for (ListObject<TYPE> *p = li.start; p != __null; p = p->next) {
      o.width(3);
      o << (++i) << ": ";
      o << (p->element) << '\n';
    }
  }
  return o;
}


////////////////////////////////////////////////////////////////////////////

template<class TYPE>
const List<TYPE> & List<TYPE>::operator= (const List<TYPE> &Li) {
  if (start != Li.start) {
    Clear(*this);
    len=Li.len;
    maxlen=len;
   
    start = copyObject<TYPE>(Li.start);  
   
    end = copyObject<TYPE>(Li.end);
    current=start;
    lastcur=__null;
   
  }
  return *this;
}

//function that copies elements, 
//support function for the copy-constructor of List

template <class TYPE> 
ListObject<TYPE> *copyObject(ListObject<TYPE> *p) {
  if (p==__null)
    return __null;
  else
    return new ListObject<TYPE>(copyObject<TYPE>(p->next), p->element);
}

////////////////////////////////////////////////////////////////////////////

#if defined(__linux) && !defined(AT_HOME)

template<class TYPE>
TYPE & First (List<TYPE> & li)
{
  li.current = li.start;
  li.lastcur = __null;
  if (li.current == __null) 
    { printf("List (First): empty list\n"); exit(1); }
  return li.current->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & Next (List<TYPE> & li)
{
  if ((li.current == __null) && (li.lastcur == __null)) return First(li);
  if (li.current == __null) return li.lastcur->element;
  li.lastcur = li.current;
  li.current = li.current->next;
  return (li.current == __null) ? li.lastcur->element : li.current->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & Last (const List<TYPE> & li) 
{
  if (li.end == __null) 
    { printf("List (Last): empty list\n"); exit(1); }
  return li.end->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & Current (const List<TYPE> & li) 
{
  if (li.current == __null) 
    { printf("List (Current): no element\n"); exit(1); }
  return li.current->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE Pop (List<TYPE> & li) 
{
  TYPE r = First(li);
  RemoveCurrent(li);
  return r;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void RemoveCurrent (List<TYPE> & li)
{
  ListObject<TYPE> *del_cur = li.current;
  if (li.current == __null) 
    { printf("List (RemoveCurrent): no element\n"); exit(1); }
  li.current = li.current->next;
  delete del_cur;
  li.len--;
  if (li.lastcur == __null) li.start = li.current;
  else li.lastcur->next = li.current;
  if (li.current == __null) li.end = li.lastcur;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void Clear (List<TYPE> &li)
{
  if( IsEmpty(li) ) return;
  while( !IsEmpty(li) )
    --li;
}




#endif // __linux

////////////////////////////////////////////////////////////////////////////

#endif // MULT_LIST_H 
