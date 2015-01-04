#ifndef DYNAMICARRAY2D_H
#define DYNAMICARRAY2D_H

template <typename T> class  DynArray2D
{
public:
    DynArray2D(int n, int m)
    {
        _n = n;
        _m = m;
    _array = new T*[n];
    for(int i = 0; i < n; i++)
       _array[i]= new T[m];

    }
    void setValue(int n, int m, T v){_array[n][m]=v;}
    T getValue(int n, int m){return _array[n][m];}
    int width(){return _n;}
    int height(){return _m;}
    ~DynArray2D(){
    for(int i=0; i < _n;i++)
        delete [] _array[i];
    delete [] _array;
    }

private:
    T **_array;
    int _n, _m;

};

#endif // DYNAMICARRAY2D_H
