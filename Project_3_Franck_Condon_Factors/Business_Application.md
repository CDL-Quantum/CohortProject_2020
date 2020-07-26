![CDL 2020 Cohort Project](../figures/CDL_logo.jpg)

# In Vino Veritas: Wine Authentication

## Step 1: Explain the technical problem you solved in this exercise
The focus of this project is to study spectroscopy from a quantum viewpoint. The Frank Condon Factor explains the intensity of vibronic transitions. Vibronic transitions are the changes in the electronic and vibrational energy of molecules when interacting with electromagnetic radiation. Spectroscopy is the study of matter and its interaction with electromagnetic radiation. Through spectroscopy and Frank Condon Factors, we can study and understand the properties of molecules. A comparison of the spectroscopic data from a sample can tell us about its temperature and chemical composition. The spectroscopic data is a fingerprint! It can be compared to other well-known sample templates to pin-down the origin and location of the sample.

To tackle the technical challenge, we studied different theoretical methods to investigate and predict chemical properties of molecules. In particular we are focusing on the so called vibronic transitions and Franck-Condon factors of molecules. First, we started by analyzing hydrogen as a proof of concept, demonstrating that our theoretical calculations line up well with experimental data. We then leveraged several sophisticated tools to expand the analysis to a variety of molecules both large and small including the Vanadium molecule. 


## Step 2: Explain or provide examples of the types of real-world problems this solution can solve

#### Astronomical Studies:
  >Astronomers leverage astronomical spectroscopy in the study of chemical compositions belonging to distant stars. The absorption lines of the sun and the emission spectra of different gases help identify the chemical constituents of the stars.
  
#### Food and Beverage:
  >The food industry currently leverages spectroscopy in the analyses required to improve the quality control of foods and beverages. Vendors can determine the falsification, adulterations, identification of origin, the origin and variety of food, beverages, and especially wine. 
  
#### Forensic Applications
  >Forensic Scientists will use spectroscopy for :
  
>Toxicology Analysis:
	Forensic analysts can identify the presence of toxic substances and their concentrations found in tissue samples and bodily fluids. Investigators can gain valuable information about the time and dosage of any poison. 
  
>Trace Evidence:
	Investigators can analyze trace elements such as fibers, glass fragments, or paint flakes. Spectroscopy can determine the exact combination of dyes used in carpet fibers and constituent components of any material analyzed. 
  
>Arson Investigations:
	investigators can analyze any residue found at a crime scene and generate accurate reports on the molecular makeup of the substance.

## Step 3: Identify at least one potential customer for this solution - ie: a business who has this problem and would consider paying to have this problem solved
The global wine market in 2018 had a valuation of 354.7 billion U.S dollars. The projected market growth is 21% by 2023 to a value of 429 billion U.S dollars.  However, the market is currently facing a wave of counterfeit products. Counter fitters are flooding the market with sub-par products, which not only takes profit away from the vineyards but also poses a threat to brand reputability. Estimates place at least 20% of the world wine supply as fake. Liv-ex a worldwide secondary market for wine with a believed approximate value of $5 billion. If as much as 20% of their wine is fake, then the counterfeit trade for their product is currently worth $1 billion. This year, with the global pandemic severely limiting wine collectors' and consumers' ability to buy wine in person directly from the source, this percentage is only expected to grow. With such large sums of money lost to fraud and the growing consumer anxiety, the wine industry is in desperate need of a solution that can verify the authenticity of their products.

The problem is that today, nearly all the anti-fraud techniques currently available are cosmetic, meaning they simply involve analyzing the visual appearance of the wine, bottle, and label. These methods are notoriously unreliable, especially when a large portion of wine counterfeiters deal in refilling authentic bottles. When it comes to analyzing the chemical composition of wine, the importance of the nearly trace amounts of compounds present vastly outweighs the major components of wine which is made up mostly of water and ethanol. So in order to accurately assesss a wine's authenticity, there needs to be a highly-sensitive and accessible solution for wine experts, investigators, collectors, and producers and this is what we would aim to accomplish. 

In our research, there are two main concerns around counterfeit wine for which we hypothesize wine industry stakeholders would pay to have a reliable solution. 
1) Collectors and restaurants pay handsomely for expensive bottles, and want to ensure beyond reasonable doubt that what they are purchasing is the real thing. 
2) An even more lucrative problem to solve is that of reputation protection. In the wine industry, it is very important that producers maintain the prestige of their brand, and a scandal involving counterfeit wine can cost them much more than just the lost revenue.

## Impact of GPL, Apache2.0 on intellectual property rights of the coder:
#### Business Application:

The codes in Project 3 (FCF_calc) use General Public License and Apache 2.0 licenses. Discuss the similarities and differences of these with respect to intellectual property rights of the coder.  What are the advantages and disadvantages of codes licensed for the public domain and those that are licensed for private use?

#### Pros and Cons of Public Domain code and Private code:
###### Note: The following is an opinion created as part of a training exercise, may not be correct and certainly does not count as a legal opinion on the matter of software Licensing. Feel free to make it an item of discussion.

#### IP and your Startup:
You have been working and hacking away at code to solve interesting problems. You may have created some of your own original code and/or also may have leveraged public domain code. Eureka! your work is now your secret sauce for your new product for your startup and you intend on protecting your IP. 

Oh, but wait ... can you actually do that? After all you have only used one GPL licensed piece of code of about 1000 lines, slightly modified. Surely your secret sauce that accounts for 100,000 lines of code, makes that modified 1000 lines insignificant, right?

Well, no, you cannot. 

If you use GPL license code in your product, your code is also under GPL. You will not be able keep your code secret. 

We do not explore the legal ramifications of various combinations here (do talk to a lawyer for your specific cases), but it is important to incorporate an IP strategy, relating to new or existing work, into your startup.

#### Licenses:
There are various types of licenses but most license applicable to open source software are designed to foster collaboration between creators with de-risking relative liability and royalty payments that could impair technological advancement. They only ask, with variations, that original creators be credited and that the new work be subject of similar constraints. Often that means one may not be able to keep their code secret. 

#### What is "Copyleft"?
Copyleft is the practice of granting the right to freely distribute and modify intellectual property with the requirement that the same rights be preserved in derivative works created from that property. 

Copyleft licenses include the **GNU General Public License (GPL)**. 

#### What are "Permissive" licenses?:
Permissive software licenses are "non-copyleft" licenses. Permissive licenses include **GNU All-permissive License, MIT License, BSD licenses, Apple Public Source License and Apache license**

#### GNU Public License (GPL):
GPL in software essentially means this: **If you modify a version of a GPL license code, your work must be under GPL also**
So in our introduction example, if your software uses any GPL code, or depends on using (not even modifying) a third-party library that is under GPL, your entire software product must be released under GPL. 

This has the following implications (among others) for your software: 
	
		First, you need to have the GPL license text with your code, indicating it is subject to GPL
    	Your source code must be made available to the public.
    	Therefore, you cannot declare your code as Intellectual Property to avoid disclosing the source

#### MIT License:
On the other hand, MIT and Apache licenses (see below for Apache) are more permissive. They require little more than attributing the original portions of the licensed code to the original developers in your own code and/or documentation. 

For example, with an MIT license (as the one used for the CDL Github, i.e. your work at CDL), 

**you can**: 
    
    re-use the code freely for your own use 
    re-use the code freely for non-commercial AND commercial re-distribution, whether in source or binary form.

**you cannot**: 
    
    Claim authorship of the software,
    Therefore, you cannot attack the original author for using or publishing his original version.

#### Apache License:

The Apache License, like MIT, is also a permissive license. However, it has stringent terms when it comes to modifications:

    You are required to explicitly list out all the modifications that you have done in the original software
    You cannot name your product in any way that hints at the product being endorsed by Apache. 

#### Advantages and Disadvantages of GPL vs MIT/Apache:

The advantages of GPL for your new original work allows you to preserve the credits to the original work and prevents others from trying to assimilate it into proprietary. It also protects you in that GPL states that no warranties are offered. 

The disadvantage as a user of original GPL work into your own work, is that you may not be able to keep your work proprietary. This is important to consider especially under a code audit for when, let's say, you are executing an exit strategy for your successful startup.

The advantage of an Apache or MIT, for the few restrictions relating to giving credit where credit is due, is that you may be able to create proprietary software whilst using third-party software. 

#### Advantages and Disadvantages of Private vs Public:
Now, regardless of the public licenses you may use, there remains the question of whether you keep your code public or private. 

Keeping your code private protects your IP at the cost of higher responsibility for your code. The buck stops with you. 
Putting you code in the public domain, if it does not hurt the business plan, can de-risk part of your venture.

It should also be noted that the rate of change in public product source code is often different that that of private products. If your code depends on third parties, you need the ability to mitigate sudden conflicting changes. Privately held software may offer a lower rate of change thus having lesser risk to your product integrity.  

There are examples where firms use both license: Private and Public. One can consider the case of MySQL AB, Red Hat and others, who offer both public and private licenses, as a means to circumvent the restrictions of GPL by offering the choice of a paid-for version or a free to use version with published source code. 

#### Conclusion

There are two major considerations with dealing with Licenses in your startup:

    What licenses will you choose for your products
    What licenses are active with the code and products you are using that have a critical impact on your business
    
Make sure your team consults the right people to make the right decisions.

    
## Step 4: Prepare a 90 second video explaining the value proposition of your innovation to this potential customer in non-technical language

[![Introduction](../figures/video.png)](https://youtu.be/vgmyN519Mew)
