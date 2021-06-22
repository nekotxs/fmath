#include <cmath>
#include <stdexcept>
#include <vector>

#ifndef FMATH_VECTOR
#define FMATH_VECTOR

namespace fmath {
    //TODO work with numeric types и в целом работу с типами, функции для 2ух и 3х мерных векторов
    template <typename T>
    class Vector {
    private:
        /**
         * @brief coords saves vector values.
         * 
         */
        std::vector<T> coords;

        void optimize();
        void optimize(const size_t &);

    public:
        Vector() = default;
        Vector(const Vector &vec);
        Vector(Vector &&vec);
        Vector(size_t size, T num);
        Vector(std::initializer_list<T> list);

        Vector &operator=(const Vector &vec);
        Vector &operator=(Vector &&vec);

        void swap(Vector &vec);
        size_t dim() const;
        T len() const;

        T &operator[](const size_t &num);
        const T operator[](const size_t &num) const;

        Vector &operator+=(const Vector &);
        Vector &operator-=(const Vector &);
        //wrong
        //Vector &operator*=(const Vector &);
        Vector &operator*=(const T &);
        Vector &operator/=(const T &);

        Vector operator+() const;
        Vector operator-() const;

        Vector operator+(Vector) const;
        Vector operator-(const Vector &) const;
        //wrong
        //Vector operator*(Vector) const;
        T operator*(const Vector &) const;

        Vector operator/(const T &) const;

        template <typename T1>
        friend Vector<T1> operator*(Vector<T1>, const T1 &);
        template <typename T1>
        friend Vector<T1> operator*(const T1 &, Vector<T1>);

        bool operator==(const Vector &) const;
        bool operator!=(const Vector &) const;

        Vector<T> vectorPow(Vector<T> second) const;
        double angle(Vector<T> second) const;

        ~Vector() = default;
    };

    template <typename T>
    Vector<T>::Vector(const Vector<T> &vec) : coords(vec.coords) {}

    template <typename T>
    Vector<T>::Vector(Vector<T> &&vec) : Vector() {
        swap(vec);
    }

    template <typename T>
    Vector<T>::Vector(size_t size, T num) : coords(size, num) {
        optimize();
    }

    template <typename T>
    Vector<T>::Vector(std::initializer_list<T> list) : coords() {
        for (auto i : list) {
            coords.push_back(i);
        }
        optimize();
    }

    template <typename T>
    Vector<T> &Vector<T>::operator=(const Vector<T> &vec) {
        coords = vec.coords;
        return *this;
    }

    template <typename T>
    Vector<T> &Vector<T>::operator=(Vector<T> &&vec) {
        swap(vec);
        return *this;
    }

    template <typename T>
    void Vector<T>::swap(Vector &vec) {
        coords.swap(vec.coords);
    }

    template <typename T>
    size_t Vector<T>::dim() const {
        return coords.size();
    }

    template <typename T>
    T Vector<T>::len() const {
        T sqsum = T();
        for (auto &i : coords) {
            sqsum += i * i;
        }
        return sqrt(sqsum);
    }

    template <typename T>
    T &Vector<T>::operator[](const size_t &num) {
        if (num < coords.size()) {
            optimize(coords.size() - num);
            return coords[num];
        }
        coords.resize(num + 1, T());
        return coords[num];
    }

    template <typename T>
    const T Vector<T>::operator[](const size_t &num) const {
        if (num < coords.size()) {
            return coords[num];
        }
        return T();
    }

    template <typename T>
    void Vector<T>::optimize() {
        while (coords.back() == T()) {
            coords.pop_back();
        }
    }

    template <typename T>
    void Vector<T>::optimize(const size_t &num) {
        for (size_t i = 0; coords.back() == T() && i != num; ++i) {
            coords.pop_back();
        }
    }

    template <typename T>
    Vector<T> &Vector<T>::operator+=(const Vector<T> &vec) {
        if (coords.size() < vec.coords.size()) {
            coords.resize(vec.coords.size(), T());
        }
        for (size_t i = 0; i != coords.size(); ++i) {
            coords[i] += vec.coords[i];
        }
        return *this;
    }

    template <typename T>
    Vector<T> &Vector<T>::operator-=(const Vector<T> &vec) {
        if (coords.size() < vec.coords.size()) {
            coords.resize(vec.coords.size(), T());
        }
        for (size_t i = 0; i != coords.size(); ++i) {
            coords[i] -= vec.coords[i];
        }
        return *this;
    }
    /*
    template <typename T>
    Vector<T> &Vector<T>::operator*=(const Vector<T> &vec) {
        if (coords.size() < vec.coords.size()) {
            coords.resize(vec.coords.size(), T());
        }
        for (size_t i = 0; i != coords.size(); ++i) {
            coords[i] *= vec.coords[i];
        }
        return *this;
    }
    */
    template <typename T>
    Vector<T> &Vector<T>::operator*=(const T &num) {
        for (auto &i : coords) {
            i *= num;
        }
        return *this;
    }

    template <typename T>
    Vector<T> &Vector<T>::operator/=(const T &num) {
        for (auto &i : coords) {
            i /= num;
        }
        return *this;
    }

    template <typename T>
    Vector<T> Vector<T>::operator+() const {
        return *this;
    }

    template <typename T>
    Vector<T> Vector<T>::operator-() const {
        Vector<T> out(*this);
        for (auto &i : out.coords) {
            i = -i;
        }
        return out;
    }

    template <typename T>
    Vector<T> Vector<T>::operator+(Vector<T> vec) const {
        vec += *this;
        return vec;
    }

    template <typename T>
    Vector<T> Vector<T>::operator-(const Vector<T> &vec) const {
        Vector<T> out(*this);
        out -= vec;
        return vec;
    }
    /*
    template <typename T>
    Vector<T> Vector<T>::operator*(Vector<T> vec) const {
        vec *= *this;
        return vec;
    }
    */
    template <typename T>
    T Vector<T>::operator*(const Vector<T> &num) const {
        T scalsum = T();
        for (size_t i = 0; i != (dim() < num.dim() ? dim() : num.dim()); ++i) {
            scalsum += coords[i] * num.coords[i];
        }
        return scalsum;
    }

    template <typename T>
    Vector<T> Vector<T>::operator/(const T &num) const {
        Vector<T> out(*this);
        for (auto &i : out.coords) {
            i /= num;
        }
        return out;
    }

    template <typename T>
    Vector<T> operator*(Vector<T> vec, const T &num) {
        for (auto &i : vec.coords) {
            i *= num;
        }
        return vec;
    }

    template <typename T>
    Vector<T> operator*(const T &num, Vector<T> vec) {
        for (auto &i : vec.coords) {
            i *= num;
        }
        return vec;
    }

    template <typename T>
    bool Vector<T>::operator==(const Vector<T> &vec) const {
        if (dim() != vec.dim())
            return false;
        for (size_t i = 0; i != dim(); ++i) {
            if (coords[i] != vec.coords[i])
                return false;
        }
        return true;
    }

    template <typename T>
    bool Vector<T>::operator!=(const Vector<T> &vec) const {
        return !(*this == vec);
    }

    template <typename T>
    Vector<T> Vector<T>::vectorPow(Vector<T> second) const {
        if (this->dim() > 3 || second.dim() > 3) {
            throw std::invalid_argument("Vector space dimension > 3");
        }
        return Vector<T>{
            (*this)[1] * second[2] - (*this)[2] * second[1],
            (*this)[0] * second[2] - (*this)[2] * second[0],
            (*this)[0] * second[1] - (*this)[1] * second[0]};
    }
    
    template <typename T>
    double Vector<T>::angle(Vector<T> second) const {
        return acos((*this * second) / (this->len() * second.len()));
    }

} // namespace fmath

#endif