import{_ as s,c as a,o as n,a7 as i}from"./chunks/framework.BUg5xOfn.js";const u=JSON.parse('{"title":"Using other functions","description":"","frontmatter":{},"headers":[],"relativePath":"article01.md","filePath":"article01.md","lastUpdated":null}'),e={name:"article01.md"},t=i(`<h1 id="Using-other-functions" tabindex="-1">Using other functions <a class="header-anchor" href="#Using-other-functions" aria-label="Permalink to &quot;Using other functions {#Using-other-functions}&quot;">​</a></h1><p>We can now calculate a transition probability matrix from a <code>LongDNA</code> sequence using the <code>transition_probability_matrix</code> and <code>initials</code> methods for a given <code>LongDNA</code> sequence:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> BioSequences, GeneFinder</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">sequence </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> dna</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;CCTCCCGGACCCTGGGCTCGGGAC&quot;</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">tpm </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> transition_probability_matrix</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sequence)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">initials </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> initials</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(sequence)</span></span>
<span class="line"></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">println</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(tpm)</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">println</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(initials)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span>4×4 Matrix{Float64}:</span></span>
<span class="line"><span> 0.0   1.0       0.0       0.0</span></span>
<span class="line"><span> 0.0   0.5       0.2       0.3</span></span>
<span class="line"><span> 0.25  0.125     0.625     0.0</span></span>
<span class="line"><span> 0.0   0.666667  0.333333  0.0</span></span>
<span class="line"><span></span></span>
<span class="line"><span>4-element Vector{Float64}:</span></span>
<span class="line"><span> 0.08695652173913043</span></span>
<span class="line"><span> 0.43478260869565216</span></span>
<span class="line"><span> 0.34782608695652173</span></span>
<span class="line"><span> 0.13043478260869565</span></span></code></pre></div><p>More conveniently, we can now use the <code>transition_model</code> method and obtain the transition probabilities and the initial distribution and build a transition model (<code>BioMarkovChain</code>):</p>`,5),p=[t];function l(o,h,c,r,d,k){return n(),a("div",null,p)}const E=s(e,[["render",l]]);export{u as __pageData,E as default};
