#include <stdio.h>
#include <stdlib.h>

#include "pointerutil.h"

void startFirstPart() {
	int firstValue;

	printf("Introduceti un numar:\n");
	scanf("%d", &firstValue);
	printf("\nDupa citire, valoarea este, dupa cum stim deja: %d.\n", firstValue);

	firstValue = 7;
	printf("I-am atribuit o alta valoare si acum este: %d.\n", firstValue);

	int *pFirst = &firstValue;
	printf("Adresa variabilei \"firstValue\" este: %p.\n", &firstValue);
	printf("O retinem intr-un pointer la int, pFirst: %p.\n", &firstValue);

	*pFirst = 5;
	printf(
		"I-am schimbat valoarea lui firstValue prin dereferentierea lui "
		"pFirst, adica modificarea valorii de la adresa %p: %d\n",
		pFirst, firstValue
	);
}

void startSecondPart() {
	int secondValue = 1 + (2 << 8) + (3 << 16) + (4 << 24);
	printf(
		"Cream o variabila, secondValue, ai carei octeti vrem sa fie de la cel"
		"mai mic, la cel mai mare: 1, 2, 3, 4.\nVa arata astfel: [4][3][2][1].\n"
	);
	printf("Valoarea decimala: %d.\n", secondValue);

	// Moduri diferite de a accesa separat cei 4 octeti
	// MODUL 1
	int leastByteMask = (1 << 8) - 1;
	printf("\nAfisarea separata a octetilor. Modul 1.\n");
	printf("---------------------------------------\n");
	printf("Cel mai mic octet din secondValue: %d.\n",
		secondValue & leastByteMask);
	printf("Al doilea octet din secondValue: %d.\n",
		(secondValue >> 8) & leastByteMask);
	printf("Al treilea octet din secondValue: %d.\n",
		(secondValue >> 16) & leastByteMask);
	printf("Cel mai mare octet din secondValue: %d.\n",
		(secondValue >> 24) & leastByteMask);

	// MODUL 2
	int mask8bits = 1 << 8;
	int mask16bits = 1 << 16;
	int mask24bits = 1 << 24;
	printf("\nAfisarea separata a octetilor. Modul 2.\n");
	printf("---------------------------------------\n");
	printf("Cel mai mic octet din secondValue: %d.\n",
		secondValue % mask8bits);
	printf("Al doilea octet din secondValue: %d.\n",
		(secondValue / mask8bits) % mask8bits);
	printf("Al treilea octet din secondValue: %d.\n",
		(secondValue / mask16bits) % mask8bits);
	printf("Cel mai mare octet din secondValue: %d.\n",
		(secondValue / mask24bits) % mask8bits);

	// MODUL 3
	char *pSecond = (char *) &secondValue;
	printf("\nAfisarea separata a octetilor. Modul 3.\n");
	printf("---------------------------------------\n");
	printf("Cel mai mic octet din secondValue: %d.\n",
		*pSecond);
	printf("Al doilea octet din secondValue: %d.\n",
		*(pSecond + 1));
	printf("Al treilea octet din secondValue: %d.\n",
		*(pSecond + 2));
	printf("Cel mai mare octet din secondValue: %d.\n",
		*(pSecond + 3));

	printf("\nNe propunem sa modificam octetul cu valoarea 3 in 71.\n");
	*(pSecond + 2) = 71;
	printf("Noua valoare decimala: %d.\n", secondValue);
}

void startThirdPart() {
	printf(
		"\nNe propunem sa citim un vector de la tastatura si sa retinem "
		"numarul de elemente pare, numarul de elemente prime si media "
		"aritmetica a acestora.\nIntroduceti lungimea vectorului:\n"
	);

	int vLength, i, v[100], evenCount, primeCount;
	float average;
	scanf("%d", &vLength);
	printf("Introduceti elementele vectorului:\n");
	for (i = 0; i < vLength; ++i) {
		scanf("%d", &v[i]);
	}

	// Functia este implementata in "pointerutil.h"
	getDataFromVector(v, vLength, &evenCount, &primeCount, &average);
	printf("Elemente pare: %d.\n", evenCount);
	printf("Elemente prime: %d.\n", primeCount);
	// Afisarea mediei se face cu 2 zecimale
	printf("Media aritmetica: %.2f.\n", average);
}

void startFourthPart() {
	printf(
		"\nUrmeaza sa intelegem notiunea de \"const\". Urmariti codul din "
		"sursa, din functia \"startFourthPart()\"."
	);

	// Keyword-ul "const" se ataseaza keyword-ului din stanga lui.
	// Daca nu exista alt keyword la stanga se ataseaza la dreapta.
	// int const* p => pointer la intreg constant
	// int *const p => pointer constant la intreg
	// const int* const p => pointer constant la intreg constant
	// const int* p => Ce e aici?

	// Ce se intampla daca incercam sa compilam codul urmator?
	// int a = 5;
	// int const* p = &a;
	// *p = 7;

	// Ce se intampla daca incercam sa compilam codul urmator?
	// int a = 5;
	// char *const p = &a;
	// *p = 7;

	//Ce se intampla daca incercam sa compilam codul urmator?
	// int a = 5;
	// char *const p = &a;
	// p ++;
	// *p = 7;
}

void startFifthPart() {
	printf("\nUrmeaza sa prezentam pointerii la functii.\n");
	// De ce am avea nevoie de pointeri la functii? Functia "qsort" din stdlib
	// primeste un pointer la o functie de comparare ca sa le sorteze in orice
	// mod isi doreste cel care o foloseste.

	// void qsort(
	// 	void *base,
	// 	size_t nitems,
	// 	size_t size,
	// 	int (*compar)(const void *, const void*)
	// )

	// Observam ca vrea sa primeasca la al 4-lea parametru un pointer la o
	// functie care returneaza un "int" si primeste ca parametri doi pointeri
	// de tip (void *), adica generici (orice), la doua variabile pe care le
	// va compara.
	// Conform http://www.cplusplus.com/reference/cstdlib/qsort/ daca functia
	// de comparare returneaza < 0, inseamna ca primul element va fi inaintea
	// celui de-al doilea in vector, la 0 sunt egale si la > 0, primul element
	// va fi pus dupa al doilea in vector 

	// Ne definim in "pointerutil.h" propria functie de comparare si invatam
	// si putin "polimorfism in C" cu ocazia acelui "void *". Functia noastra
	// va sorta DESCRESCATOR, deci daca la < 0 ele vor ramane in pozitia
	// initiala (primul inainte de al doilea) atunci vom returna diferenta
	// dintre valoarea celui de-al doilea element si a primului.
	// primul = 9, al doilea = 7. 7 - 9 = -2 < 0 deci vor ramane in ordinea
	// (9, 7), adica DESCRESCATOR
	
	int v[] = {1, 7, 29, 2, 12, 4, 0};
	printf("Inainte de sortare:\n");
	printVector(v, 7);
	
	qsort(v, 7, sizeof(v[0]) /* 4 */, descendingComparingFunction);
	printf("Dupa sortare:\n");
	printVector(v, 7);

}

int main() {
	startFirstPart();
	startSecondPart();
	startThirdPart();
	startFourthPart();
	startFifthPart();
	return 0;
}