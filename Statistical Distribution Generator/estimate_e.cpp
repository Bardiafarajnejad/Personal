//
//  
//  estimate_e
//
//  Created by Bardia Farajnejad on 7/17/20.
//  Copyright © 2020 Bardia Farajnejad. All rights reserved.
//

#include <iostream>
#include <iomanip>
using namespace std;

int main()
{
    // insert code here...
    int number1; // first integer read from user
    double sum=1; //   I AM MANUALLY ADDING +1 THE FIRST TERM HERE && sum must be a double
    double term;
    cout << "Enter the number of terms in the sum (highter accuracy--> more terms): ";
    cin >> number1; // read values from user

    for (int i=2; i<=number1; i=i+1)
    {
    int factorial=1;
    int count=1;
        while (count < i)
        {
            factorial=factorial * (count);
            count=count+1;
        }
        term = static_cast< double >(1)/factorial;
        sum=sum+term;
        cout << "\nthe term is " << setprecision(100) << fixed << term;
        cout << "\nthe sum is " << setprecision(100) << fixed << sum;
    }
    cout << "\nThe sum of the first " << number1 << " terms is " << setprecision(30) << fixed << sum << endl;
} // end main

