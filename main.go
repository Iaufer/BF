package main

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
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

func Weight_int(a uint32) uint32 {
	x := a
	x = x - ((x >> 1) & 0x55555555) // ∧ 
	x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
	x = (x + (x >> 4)) & 0x0F0F0F0F
	x = x + (x >> 8)
	x = x + (x >> 16)

	return x & 0x3F
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
	count := 0
	for i := 0; i < bf.nw; i++ {
		if bf.f[i] == 0 {
			count++
		}
	}

	if count == bf.nw {
		fmt.Println(0)
		return
	}

	f := 0
	// fmt.Println(bf.n)
	if bf.f[0]&uint32(1) == 1 {
		print("1")
		f = 0
	} else {
		print("0")
		f = 1
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
						print("x", j)
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

func WH_tr(bf *BF) []int {
	//убрать умножение
	cf := make([]int, 1<<bf.n) //128 секунд для 28 перменных, для 27 20 секн
	// cf := make([]int, 0) // аллокации и append

	var mask, size_j uint32 = 1, 0

	if bf.n > 5 { // для того чтобы чтобы определить кочичество бит в ячейки если количество перменных меньше  или равно 5(4. 8. 16. 32), а если перменных больше 5 то бит в ячейки всегда 32
		size_j = 32
	} else {
		size_j = 1 << bf.n
	}

	for i := uint32(0); i < uint32(bf.nw); i++ {
		for j := uint32(0); j < size_j; j++ {
			if bf.f[i]&mask == 0 {
				// cf = append(cf, 1)
				cf[j+i<<5] = 1
			} else {
				cf[j+i<<5] = -1
				// cf = append(cf, -1)
			}
			mask = mask << 1
		}
		mask = 1
		// count++
	}

	step := 1

	size_cf := len(cf)

	for step < size_cf {
		for i := 0; i < size_cf; i += step << 1 {
			for j := i; j < i+step; j++ {
				a := cf[j]
				b := cf[j+step]
				cf[j] = a + b
				cf[j+step] = a - b
			}
		}
		step = step << 1
		// step *= 2
	}

	// fmt.Println(cf, len(cf))

	return cf
}

//11000000 -
//найти все Линейный функцииБ фиктивные
func (bf *BF) cor(cf []int) int {

	for i := uint32(1); i < (1 << bf.n); i++ {
		next, end := next_comb(i, bf, cf)
		if next == -1 {
			return -1
		}

		if end == -1 {
			return next
		}

	}
	return int(bf.n)

}

func next_comb(_i uint32, bf *BF, cf []int) (int, int) {
	var a, b, c, k, n uint32 = 0, 0, 0, _i, bf.n

	a = ((1 << k) - 1) << (n - k)

	if uint32(cf[int(a)]) != 0 {
		// println("cf ", (cf[int(a)]), a, k)
		return int(k - 1), -1
	}

	if k == bf.n { // чтобы вовремя остановаиться
		// fmt.Println("stop", Weight_int(a), bf.n)
		return int(k), -1
	}

	for i := uint32(0); i < bf.n-1; i++ {
		b = (a + 1) & a
		c = Weight_int((b-1)^a) - 2
		if uint32(cf[int(a)]) != 0 {
			// fmt.Println(Weight_int(a)-1, cf[a])
			return int(k - 1), -1
		}
		a = (((((a + 1) ^ a) << 1) + 1) << c) ^ b
		// println(a)

	}

	return 1, 0
}

func nonBF(bf *BF, _cf []int) int {
	cf := make([]int, len(_cf))
	copy(cf, _cf)
	// for i := 0; i < len(cf); i++ {
	// 	if cf[i] < 0 {
	// 		cf[i] *= -1
	// 	}
	// }
	// // fmt.Println(cf)

	// sort.Ints(cf)

	maxAbs := cf[0]
	for i := 1; i < len(cf); i++ {
	    if cf[i] < 0 {
	        cf[i] *= -1
	    }
	    if cf[i] > maxAbs {
	        maxAbs = cf[i]
	    }
	}

	// fmt.Println("cf: ", cf[len(cf)-1])
	return (1 << (bf.n - 1)) - (maxAbs/2)
	
	// return (1 << (bf.n - 1)) - (cf[len(cf)-1])/2
}

////Best Affine Approximations
func BAA(bf *BF, cf []int) {
	fmt.Println("Наилучшее аффиное приближение: ")
	for i := uint32(0); i < (1 << bf.n); i++ {
		if cf[i] < 0 {
			fmt.Print("1")
			print_anf_BAA(i, bf.n)
		} else if cf[i] > 0 {
			fmt.Print("0")
			print_anf_BAA(i, bf.n)
		}
	}
}

func print_anf_BAA(_i, n uint32) {
	var mask uint32 = 1
	for i := uint32(0); i < n; i++ {
		if _i&mask != 0 {
			fmt.Print(" + x", i)
		}
		mask = mask << 1
	}
	fmt.Println()
}

//∆f (x)
func auto_cor(cf []int) []int {
	_cf := make([]int, len(cf))
	for i := 0; i < len(cf); i++ {
		_cf[i] = cf[i] * cf[i]
	}


	step := 1

	for step < len(_cf) {
		for i := 0; i < len(_cf); i += step * 2 {
			for j := i; j < i+step; j++ {
				a := _cf[j]
				b := _cf[j+step]
				_cf[j] = a + b
				_cf[j+step] = a - b
			}
		}
		step *= 2
	}

	fmt.Println(_cf)

	d := uint32(math.Log2(float64(len(_cf)) + 1))

	for i := 0; i < len(_cf); i++ {
		_cf[i] = _cf[i] / (1 << d)
	}

	return _cf
}

//фиктивные переменные(dummy variables)
func d_var(_bf *BF) {
	bf := constr_Copy(_bf)
	meb := Mebius(bf)

	// print_BF1(meb)
	var mask, tmp uint32 = 1, 0

	for i := uint32(0); i < (1 << meb.n); i++ {
		if meb.f[i/32]&mask != 0 {
			tmp |= i
			// println(mask)
		}
		mask = mask << 1
	}

	if tmp == (1<<bf.n)-1 { //линейносchr
		fmt.Println("Фиктивных перeменных нет!")
		return
	}

	mask = 1

	for i := int(meb.n) - 1; i >= 0; i-- {
		if tmp&mask == 0 {
			// fmt.Sprintf("x %d", )
			fmt.Print("x", i)
		}
		mask = mask << 1
	}
	println()

	// for i := uint32(0); i < meb.n; i++ {
	// 	if tmp&mask == 0 {
	// 		fmt.Println("x", i)
	// 	}
	// 	mask = mask << 1
	// }

	// println(tmp)
}

//Линейный переменные
func lin_var(_bf *BF) {
	bf := constr_Copy(_bf)
	meb := Mebius(bf)

	// print_BF1(meb)
	var mask, tmp uint32 = 1, 0

	for i := uint32(0); i < (1 << meb.n); i++ {
		if meb.f[i/32]&(mask) != 0 && Weight_int(i) == 1 {
			tmp |= i
		}
		mask = mask << 1 //(2<<i / вместо маски
	}

	mask = 1

	_tmp := uint32(0)

	for i := uint32(0); i < (1 << meb.n); i++ {
		if meb.f[i/32]&(mask) != 0 && Weight_int(i) != 1 {
			_tmp |= i
		}
		mask = mask << 1 //(2<<i / вместо маски
	}
	mask = 1

	for i := uint32(0); i < meb.n; i++ {
		_tmp ^= mask
		mask <<= 1
	}

	tmp = tmp & _tmp
	mask = 1

	// println(tmp)
	// if tmp == 0 {
	// 	fmt.Println("Все переменные нелинейны!")
	// 	return
	// }

	for i := int(meb.n) - 1; i >= 0; i-- {
		if tmp&mask != 0 {
			fmt.Printf("x%d\n", i)
		}
		mask = mask << 1
	}
	println()

}

func main() {
	rand.Seed(time.Now().UnixNano())
	a := "11111111"
	bf := BF_constructorch(a)
	// bf := BF_constructor(7, 2)
	// //3
	wh_tr := WH_tr(bf)
	fmt.Println("Преобразование Уолша – Адамара: ", wh_tr)
	meb := Mebius(bf)
	fmt.Println("Mebius:")
	print_BF1(meb)

	cor := bf.cor(WH_tr(bf))
	fmt.Println("Максимальный порядок корреляционной иммунности функции: ", cor)

	// // //4
	nonbf := nonBF(bf, wh_tr)
	fmt.Println("Нелинейность функции: ", nonbf)

	BAA(bf, wh_tr) //NAP

	fmt.Println("Фиктивные переменные: ")
	d_var(bf)

	// fmt.Println("Линейныe переменные: ")
	// lin_var(bf)

}
