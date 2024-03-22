package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"
)

type BF struct {
	f    []uint32
	n    uint32 //количество переменных
	nw   int    //число ячеек
	nbit int    // кол-во бит
}

func BF_constructor(n, typ uint32) *BF {
	var bf BF

	switch typ {
	case 0:
		bf.n = n
		bf.nw = ((1 << n) + 31) >> 5
		bf.f = make([]uint32, bf.nw)
	case 1:
		bf.n = n
		bf.nw = ((1 << n) + 31) >> 5
		bf.f = make([]uint32, bf.nw)
		if bf.n >= 5 {
			for i := 0; i < bf.nw; i++ {
				bf.f[i] = ^uint32(0)
			}
		} else {
			bf.f[0] = (1<<(uint32(1<<bf.n)) - 1)
		}
	case 2:

		bf.n = n
		bf.nw = ((1 << n) + 31) >> 5
		bf.f = make([]uint32, bf.nw)

		if bf.n >= 5 {
			for i := 0; i < bf.nw; i++ {
				bf.f[i] = rand.Uint32()
			}
		} else {
			bf.f[0] = rand.Uint32() >> (uint32(32) - uint32(1<<bf.n)) //в

		}
	}

	return &bf
}

func BF_constructorch(vec string) *BF {
	var bf BF

	if len(vec) != 0 && (len(vec)&(len(vec)-1)) == 0 {

		bf.nbit = len(vec)

		bf.nw = ((bf.nbit - 1) / 32) + 1

		bf.f = make([]uint32, bf.nw)

		bf.n = uint32(math.Log2(float64(bf.nbit) + 1))

		for i := 0; i < bf.nbit; i++ {
			if vec[i] != '0' {
				bf.f[i/32] = bf.f[i/32] | (1 << (i % 32))
			}
		}

		return &bf
	} else {
		println("Invalid")
		bf.f = nil
		bf.n = 0
		bf.nw = 0

		return &bf
	}
}

func print_BF1(bf *BF) {
	size := int(math.Pow(2, float64(bf.n)))

	for i := 0; i < size; i++ {
		if (bf.f[i/32] & (1 << (i % 32))) == 0 {
			print(0)
		} else {
			print(1)
		}

	}
	println()
}

func (bf *BF) Equal(other *BF) bool {
	if bf.nbit != other.nbit {
		return false
	}
	for i := 0; i < bf.nw; i++ {
		if bf.f[i] != other.f[i] {
			return false
		}
	}
	return true
}

func (bf *BF) Giving(other *BF) *BF {

	if bf != other {
		bf.nw = other.nw
		bf.nbit = other.nbit

		bf.f = nil //освобождение памяти

		bf.f = make([]uint32, bf.nw)

		for i := 0; i < bf.nw; i++ {
			bf.f[i] = other.f[i]
		}
	}

	return bf
}

func Weight(bf *BF) uint32 {

	var sum uint32 = 0
	for i := 0; i < len(bf.f); i++ {
		var x uint32 = bf.f[i]

		x = x - ((x >> 1) & 0x55555555) // ∧ - конънк
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
		x = (x + (x >> 4)) & 0x0F0F0F0F
		x = x + (x >> 8)
		x = x + (x >> 16)

		sum += (x & 0x3F)
	}

	return sum
}

func Kn(bf *BF) float64 { //знаечние Kn

	wei := Weight(bf)
	n := 1 << bf.n

	return float64(wei) / float64(n)
}

func constr_Copy(bf *BF) *BF {
	newBF := &BF{}

	newBF.n = bf.n
	newBF.nbit = bf.nbit
	newBF.nw = bf.nw

	newBF.f = make([]uint32, bf.nw)

	for i := 0; i < bf.nw; i++ {
		newBF.f[i] = bf.f[i]
	}

	return newBF
}

func Mebius(_bf *BF) *BF {
	bf := constr_Copy(_bf)

	for i := 0; i < bf.nw; i++ {
		g := bf.f[i]
		g = g ^ ((g << 1) & 0xaaaaaaaa)
		g = g ^ ((g << 2) & 0xcccccccc)
		g = g ^ ((g << 4) & 0xf0f0f0f0)
		g = g ^ ((g << 8) & 0xff00ff00)
		g = g ^ ((g << 16) & 0xffff0000)
		bf.f[i] = g
	}
	step := 1

	for step < bf.nw {
		for i := 0; i < bf.nw; i += step * 2 {
			for j := i; j < i+step; j++ {
				bf.f[j+step] = bf.f[j+step] ^ bf.f[j]
			}
		}
		step *= 2
	}
	return bf
}

func ANF_print(bf *BF) {
	f := 0
	fmt.Println(bf.n)
	if bf.f[0]&uint32(1) == 1 {
		print("1")
		f = 0
	}
	g := bf.f[0]

	ch := uint32(0) //для того чтобы обновлять маску для каждой новой ячейки

	var mask, _mask1 uint32 = 2, 1
	for i := uint32(1); i < (1 << bf.n); i++ {
		g = bf.f[i/32]
		if ch != g && i != 1 {
			// fmt.Println("Обновить маску!", i, 1<<bf.n) i != 1 нужно из-за самой первой итерации
			mask = 1

		}
		ch = g
		_mask1 = 1

		for j := 0; j < int(bf.n); j++ {
			if g&mask != 0 {
				if i&_mask1 != 0 {
					if f == 0 {
						print(" + ")
						f = 1
					}
					if f == 1 {
						print("x", j+1)
					}
				}
			}
			_mask1 = _mask1 << 1
		}

		// j := int(bf.n) - 1

		// for j >= 0 {
		// 	if g&mask != 0 {
		// 		if i&_mask1 != 0 {
		// 			if f == 0 {
		// 				print(" + ")
		// 				f = 1
		// 			}
		// 			if f == 1 {
		// 				print("x", j)
		// 			}
		// 		}
		// 	}
		// 	j--
		// 	_mask1 = _mask1 << 1
		// }
		f = 0
		mask = mask << 1
	}
	println()
}

func getDeegre(bf *BF) int {
	g := bf.f[0]

	ch := uint32(0) //для того чтобы обновлять маску для каждой новой ячейки

	tmp := 0

	var mask, _mask1 uint32 = 2, 1
	count := 0
	for i := uint32(1); i < (1 << bf.n); i++ {
		count = 0
		g = bf.f[i/32]
		if ch != g && i != 1 {
			// fmt.Println("Обновить маску!", i, 1<<bf.n) i != 1 нужно из-за самой первой итерации
			mask = 1

		}
		ch = g
		_mask1 = 1

		for j := 0; j < int(bf.n); j++ {
			if g&mask != 0 {
				if i&_mask1 != 0 {
					count++
				}
			}
			_mask1 = _mask1 << 1
		}
		mask = mask << 1

		if count > tmp {
			tmp = count
		}
	}

	// println("Степень: ", tmp)

	return tmp
}

func main() {
	rand.Seed(time.Now().UnixNano())
	// a := "1010110110101101101011011010110110101101101011011010110110101101"
	// a := "1011010101101100"
	// bf := BF_constructorch(a)
	bf := BF_constructor(4, 2)

	start := time.Now()

	meb := Mebius(bf)
	// meb := Mebius(bf)

	end := time.Now()

	fmt.Println("Время работы преобразование Мебиуса: ", end.Sub(start).Milliseconds())
	d := getDeegre(meb)

	fmt.Println("Степень функции: ", d)

	// // print_BF1(bf)
	// // print_BF1(meb)

	if bf.Equal(Mebius(meb)) {
		println("Преобразование мебиуса работает корректно!")
	} else {
		println("Не работает преобразование")
	}
}
