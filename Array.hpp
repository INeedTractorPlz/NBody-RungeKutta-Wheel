#include <cstddef>
#include<iostream>

#include <boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include<boost/shared_ptr.hpp>

using namespace boost::numeric::ublas;
using namespace boost;

/*В данном hpp файле реализован класс массивов Array со следующими особенностями.
1)Array<T> A(B); Происходит копирование Умного Указателя на память, хранящую элементы B в A.
Т.е. изменение элементов А изменит элементы B и наоборот.
Это полезно при возвращении больших массивов из функции, если массив создан внутри функции.
Т.к. копирование указателя, а не элементов приведёт к серьезной экономии времени, а умный указатель
предотвратит удаление памяти, на которую указывал массив, созданный в функции.
При использовании оператора копирования следует также помнить, что память выделенная под массив затрётся
только после окончания времени существования последней копии этого массива(не поэлементной копии).
2)Array<T> A(B.copy()); Копирование поэлементное.
3)A=B; Также поэлементное копирование и затирание памяти из A(при условии, что это последняя копия).
4)A.subjoin(B); Обычное копирование(умного указателя) массива B в следующую часть массива A.
4.1)Array<Array<T> > A;
  Array<T> B;
    -//- Массивы A и B как-то заполнены.
    B.subjoin(A[0]); Не поэлементное копирование.
    В таких более сложных случаях можно забыть, что копирование не поэлементное и массив B остаётся связан
с массивом A[0].
5)A=Array<T>(B); Здесь копирование поэлементное, т.к. используется оператор присваивания.
    Если при задании переменной A не сделать её непоэлементной копией массива B,
    то потом сделать это с помощью конструктора копирования уже нельзя.
  Для этого используется:
  A.mirror(B);
*/
template<typename T>
struct Array;

template<typename U>
struct Arr
{
    // Список операций:
    //
    explicit Arr(size_t size, const U& value = U())
    //   конструктор класса, который создает
    //   Arr размера size, заполненный значениями
    //   value типа U. Считаем, что у типа U есть
    //   конструктор, который можно вызвать без
    //   без параметров, либо он ему не нужен.
    : size_(size)    {
        data_= static_cast<U*>(operator new[] (size * sizeof(U)));
        for(unsigned i=0; i< size_; i++)
            new (data_+i)U(value);
    }
    Arr(const Arr<U> &A)
    //   Конструктор копирования, который копирует поэлементно первую часть массива А и все последующие,
    //   если они есть, в первую и последующие части нашего массива соответственно.
    :    size_(A.size_) {
        data_= static_cast<U*>(operator new[] (size_ * sizeof(U)));
        for(unsigned i=0; i< size_; i++)
            new (data_+i)U(A.data_[i]);
        if(A.data_next!=0)
            data_next.reset(new Array<U>((A.data_next)->copy()));
    }
    Arr(){
        data_=0;
        size_=0;
    }   
    ~Arr(){    
    //   деструктор
        for (int i = 0; i < size_; i++)
            data_[i].~U();
        operator delete [] (data_);    
    }    
    Arr& operator=(const Arr<U> &A){
    //   Оператор присваивания.
    //  Удаляем содержимое нашего массива, выделяем новую память и заполняем значениями из массива А.
    //  Тоже самое делаем для следующей части массива, используя конструктур копирования
    //   и функцию copy, определённые в классе Array. 
        if(this->data_!=A.data_){
            this->~Arr();
            this->size_=A.size_;
            data_= static_cast<U*>(operator new[] (size_ * sizeof(U)));
            for(unsigned i=0; i< size_; i++)
                new (data_+i)U(A.data_[i]);     
        }
        if(data_next!=0)
            *data_next=*(A.data_next);
        else
            if(A.data_next!=0)
                data_next.reset(new Array<U>((A.data_next)->copy()));
        return *this;
    }    
    void subjoin(size_t size, const U& value = U()){
    //  Добавляем новую часть массива размера size, заполненную значениями value.
        if(data_next==0)
            data_next.reset(new Array<U>(size, value));
        else
            data_next->subjoin(size, value);
    }
    void subjoin(const Array<U>& A){
    //  Добавляем новую часть массива и копируем в неё Array<U> A,
    //  стоит обратить внимание на то, что копирование происходит без использования функции copy,
    // а значит копируется лишь указатель, и при изменении элементов массива A будут меняться
    // и элементы нашего массива, и наоборот. Чтобы избежать подобного необходимо использовать
    //  функцию copy при передаче аргумента.
        data_next.reset(new Array<U>(A)); 
    }
    size_t size() const{
    //  Возвращает полный размер нашего массива, т.е. сумму размеров всех его частей.
        if(data_next==0)    
            return this->size_;
        else
            return size_+data_next->size();
    }
    U& operator[](size_t i){
    //  Возвращает элемент массива по номеру, игнорируя разбиение на части.
        if(i<size_)
            return data_[i];
        else
            return (*data_next)[i-size_];
    }
    const U& operator[](size_t i) const{
    //  Возвращает элемент массива по номеру, игнорируя разбиение на части.
        if(i<size_)
            return data_[i];
        else
            return (*data_next)[i-size_];
    }
    
    U *data_;   //  Указатель на память, в которой хранится первая часть данного массива.
    shared_ptr<Array<U> > data_next;    // Умный указатель на массив, являющейся второй частью данного.
    size_t size_;   // Размер первой части данного массива.
};


// Обёртка над Arr, позволяющая копировать только указатели на память, а не сами элементы массива,
//  без явного использования указателей.
template<typename T>
struct Array{
    T& operator[](size_t i){return (*array)[i];}
    const T& operator[](size_t i) const{return (*array)[i];}
    size_t size() const{return (this->array)->size();}
    explicit Array(size_t size, const T& value = T()) : array(new Arr<T>(size,value)){}
    Array(Arr<T>* copy) : array(copy){ } // Конструктор класса, получающий на вход обычный указатель на Arr,
    // и конвертирующий его в умный, записывая в array. Очевидно, изменение элементов copy приведет к
    // изменению элементов текущего массива и наоборот.
    Array(const Array<T>& copy) : array(copy.array){ } //Конструктор, копирующий умный указатель на Arr
    // в array. Очевидно, изменение элементов copy приведет к
    // изменению элементов текущего массива и наоборот.    
    Array() : array(shared_ptr<Arr<T> >(new Arr<T>())) {} //Конструктор, создающий указатель на
    // пустой Arr.

    Array& operator=(const Array<T> &A){
    // Оператор присваивания.
    // Копирование массива A в текущий с использованием констуктора Arr, т.е. поэлементно.
    // Освобождать память не нужно, умный указатель справится с этим сам.
        this->array=shared_ptr<Arr<T> >(new Arr<T>(*(A.array)));
        return *this;
    }
    Array<T> copy(){
    // Копирование массива A в возвращаемый с использованием констуктора Arr, т.е. поэлементно.
    // Освобождать память не нужно, умный указатель справится с этим сам.
        return Array<T>(new Arr<T>(*(this->array))); 
    }
    void mirror(const Array<T>& A){
    //  Копирование указателей. Очевидно, изменение элементов А приведет к
    // изменению элементов текущего массива и наоборот.
        this->array=A.array;
    }
    void subjoin(size_t size, const T& value = T()){
        array->subjoin(size, value);
    }
    void subjoin(const Array<T>& A){
        (array->data_next).reset(new Array<T>(A)); 
    }

    shared_ptr<Arr<T> > array;  //  Умный указатель на массив Arr.
};

//  В случае массива массива конструктор класса работает некорректно. А именно создаёт массив
//  указателей на один и тот же массив. Проблема в том, что подавая в value массив Array,
//  мы будем использовать конструктор копирования Array, который копирует только указатель, а не элементы.

template<typename T>
struct Array<Array<T> >{
    explicit Array(size_t size, Array<T> value) : array(new Arr<Array<T> >(size,value)){
        for(unsigned i=1;i<size;++i)
            array->data_[i]=value.copy();// Поэлементно копируем массив value в элементы нашего массива.
    }
//  Остальное всё тоже самое. Скопировать пришлось из-за особенностей частичной специализации.
    Array(Arr<Array<T> >* copy) : array(copy){ }
    Array(const Array<Array<T> >& copy) : array(copy.array){ }    
    Array() : array(shared_ptr<Arr<Array<T> > >(new Arr<Array<T> >)) {}

    Array& operator=(const Array<Array<T> > &A){
        this->array=shared_ptr<Arr<Array<T> > >(new Arr<Array<T> >(*(A.array)));
        return *this;
    }

    Array<Array<T> > copy(){
        return Array<Array<T> >(new Arr<Array<T> >(*(this->array))); 
    }


    void mirror(const Array<Array<T> >& A){
        this->array=A.array;
    }

    void subjoin(const Array<Array<T> >& A){
        (array->data_next).reset(new Array<Array<T> >(A)); 
    }

    void subjoin(size_t size, const Array<T>& value){
        if(array->data_next==0)
            (array->data_next).reset(new Array<Array<T> >(size, value));
        else
            (array->data_next)->subjoin(size, value);
    }

    Array<T>& operator[](size_t i){return (*array)[i];}
    const Array<T>& operator[](size_t i) const{return (*array)[i];}
    size_t size() const{return (this->array)->size();}
    
    shared_ptr<Arr<Array<T> > > array;

};

//  Дальше идёт реализация различных операций с массивами.

template<typename T>
Array<T> operator+(Array<T> left, Array<T> right){
    if(left.size()!=right.size()){
        std::cout << "Разные размеры складываемых массивов!!!" << std::endl;
        exit(0);
    }
    Array<T> result(left.copy());
    for(unsigned i=0;i<result.size();++i)
        result[i]+=right[i];
    return result;
}


template<typename T>
Array<T> operator+=(Array<T> left, Array<T> right){
    left=left+right;
    return left;
}


template<typename T, typename U>
Array<T> operator*(Array<T> left, const U &right){
    Array<T> result(left.copy());
    for(unsigned i=0;i<left.size();++i)
        result[i]*=right;
    return result;
}


template<typename T, typename U>
Array<T> operator*(const U &right, Array<T> left){
    Array<T> result(left.copy());
    for(unsigned i=0;i<left.size();++i)
        result[i]*=right;
    return result;
}

template<typename T, typename U>
Array<T> operator/(Array<T> left, const U &right){
    Array<T> result(left.copy());
    for(unsigned i=0;i<left.size();++i)
        result[i]/=right;
    return result;
}


template<typename T, typename U>
Array<T> operator/(const U &right, Array<T> left){
    Array<T> result(left.copy());
    for(unsigned i=0;i<left.size();++i)
        result[i]/=right;
    return result;
}

template<typename T, typename U>
Array<T> operator*=(Array<T> left, const U &right){
    left=left*right;
    return left;
}
template<typename T, typename U>
Array<T> operator/=(Array<T> left, const U &right){
    left=left/right;
    return left;
}