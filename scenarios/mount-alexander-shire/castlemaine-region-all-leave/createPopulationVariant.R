gz<-'../castlemaine-region/population-archetypes.csv.gz'
con<-gzfile(gz,'rt')
orig<-read.csv(con,header=T,sep=',',stringsAsFactors = F,strip.white = T)
close(con)

df<-orig

# Simplofy by removing any columns that are inconsequentiual (not used by the model)
df$Id<-NULL
df$AgentId<-NULL
df$Age<-NULL
df$Archetypes.Age<-NULL
df$Archetypes.Household<-NULL
df$EZI_ADD<-NULL
df$Gender<-NULL
df$HouseholdId<-NULL
df$PrimaryFamilyType<-NULL
df$SA1_7DIGCODE<-NULL
df$SA2_MAINCODE<-NULL

# Adjust some other values
df$HasDependents<-"false"
df$WillGoHomeAfterVisitingDependents<-"false"
df$WillGoHomeBeforeLeaving<-"false"
df$WillStay<-"false"
df$ResponseThresholdInitial<-0
df$ResponseThresholdFinal<-0
df$ImpactFromFireDangerIndexRating<-0
df$ImpactFromMessageAdvice<-0
df$ImpactFromMessageEmergencyWarning<-0
df$ImpactFromMessageEvacuateNow<-1
df$ImpactFromMessageRespondersAttending<-0
df$ImpactFromMessageWatchAndAct<-0
df$ImpactFromImmersionInSmoke<-0
df$ImpactFromVisibleEmbers<-0
df$ImpactFromVisibleFire<-0
df$ImpactFromVisibleResponders<-0
df$ImpactFromVisibleSmoke<-0
df$ImpactFromSocialMessage<-0
df$LagTimeInMinsForInitialResponse<-1
df$LagTimeInMinsForFinalResponse<-1
df$InvacLocationPreference<-df$EvacLocationPreference
df$HasDependentsAtLocation<-""

# Write out
con <- gzfile('./population-archetypes.csv.gz')
write.csv(df, con, row.names=FALSE, quote=TRUE)
